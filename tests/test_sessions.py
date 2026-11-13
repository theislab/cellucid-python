"""Session format 4 in Python: the golden history read exactly, and applied.

The shared golden fixture (``tests/fixtures/history_golden``, written by the
web repository's ``scripts/generate-history-golden.mjs`` with the history
model's own constructors; its home is the web repository's
``tests/fixtures/``, and cellucid-python and cellucid-r each hold a copy the
tests compare with it byte for byte) holds every kind of record the model has: trackers (a replay,
one removed), branches (a nested fork, a locked and an archived one), turns
(a notebook's, an agent's with its summary), notes, a bookmark, a paused
tracker, squashed steps (a paused period; a squash in mid-branch carrying a
note and a bookmark, with a fork after it), a deleted step and its journal
entry, a flagged step, a stale result. ``history_golden.expected.json`` says what reading it must give, and
the apply plans of a sequence of applies on one object: these tests compare
field by field. Synthesized files (``_session_files.py``) test every other
rule: what the share left out, damaged and newer files, and each apply rule
on its own.
"""

from __future__ import annotations

import base64
import copy
import hashlib
import json
import math
import re
import struct
from pathlib import Path
from typing import Any

import _session_files as files
import anndata as ad
import numpy as np
import pandas as pd
import pytest
from _web_root import assert_shared_file, web_file

import cellucid
from cellucid import (
    ApplyPlan,
    ApplyReport,
    CellucidSession,
    apply_cellucid_session_to_anndata,
    read_session,
)
from cellucid.anndata_session._apply import (
    GENE_SETS_UNS_ONLY_NOTE,
    ORIGINAL_POSITION_MESSAGE,
    PRESENTATION_NOTE,
    cell_order_digest,
    gene_set_columns,
    history_applied,
)
from cellucid.anndata_session._colors import COLORS_NOTE
from cellucid.history import NOTE_COLUMNS, STEP_COLUMNS, Diff
from cellucid.session import (
    SESSION_SIGNATURE,
    _hash_value,
    _js_number,
    _js_string,
    _js_stringify,
    newer_session_message,
    package_version,
)

FIXTURES = Path(__file__).parent / "fixtures" / "history_golden"
GOLDEN = FIXTURES / "history_golden.cellucid-session"
EXPECTED_JSON = FIXTURES / "history_golden.expected.json"
#: The golden's current state only, as the web app's Share › Current state only writes it.
GOLDEN_STATE = FIXTURES / "history_golden_state.cellucid-session"
#: The golden session's SHA-256, which names the file in a session, a plan and an apply record.
GOLDEN_SHA256 = hashlib.sha256(GOLDEN.read_bytes()).hexdigest()
EXPECTED = json.loads(EXPECTED_JSON.read_text(encoding="utf-8"))
T1, T2, T3 = (tracker["id"] for tracker in EXPECTED["trackers"])


@pytest.fixture(scope="module")
def golden() -> CellucidSession:
    return read_session(GOLDEN)


def golden_adata() -> ad.AnnData:
    """The dataset the golden session was saved on, as expected.json describes it."""
    dataset = EXPECTED["dataset"]
    obs = pd.DataFrame(
        {
            name: pd.Categorical.from_codes(column["codes"], categories=column["categories"])
            for name, column in dataset["obs"].items()
        },
        index=pd.Index(dataset["obs_names"]),
    )
    adata = ad.AnnData(
        X=np.zeros((dataset["n_obs"], len(dataset["var_names"])), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(dataset["var_names"])),
    )
    adata.obsm["X_umap"] = np.asarray(dataset["obsm"]["X_umap"], dtype=np.float32)
    return adata


def _rows(plan: ApplyPlan) -> list[dict[str, Any]]:
    return [
        {"kind": row.kind, "domain": row.domain, "name": row.name, "column": row.column, "from": list(row.from_),
         "cells": row.cells, "status": row.status, "note": row.note}
        for row in plan.rows
    ]


# =============================================================================
# The golden fixture, read
# =============================================================================


@pytest.mark.parametrize(
    ("ours", "theirs"),
    [
        ("history_golden.cellucid-session", "history-golden.cellucid-session"),
        ("history_golden.expected.json", "history-golden.expected.json"),
        ("history_golden_state.cellucid-session", "history-golden-state.cellucid-session"),
    ],
)
def test_the_golden_files_are_the_web_apps_byte_for_byte(ours: str, theirs: str) -> None:
    assert_shared_file(f"tests/fixtures/history_golden/{ours}", f"tests/fixtures/{theirs}")


def test_the_expected_reading_names_the_golden_session() -> None:
    assert EXPECTED["sha256"] == GOLDEN_SHA256


def test_the_web_apps_current_state_only_reads_as_one_state(golden: CellucidSession) -> None:
    """Share › Current state only: one tracker of the opening and one squashed step, nobody named, no time, no note."""
    state = read_session(GOLDEN_STATE)
    assert (state.scope.content, state.scope.times, state.scope.notes, state.scope.chats) == ("state", "none", False, False)
    assert state.created_at == "1970-01-01T00:00:00.000Z"
    (tracker,) = state.trackers
    assert tracker.created_at is None
    root, shared = tracker.step_records
    assert (root.kind, root.by.origin, shared.kind, shared.squashed.cause, shared.by.origin) == (
        "root", "startup", "squashed", "share", "user"
    )
    assert tracker.head_id == shared.id
    # What the golden shows, value for value; a user field's creation time is left out as every time is.
    shown = tracker.step(shared.id).state()
    was = golden.tracker().step(golden.tracker().head_id).state()
    assert set(shown) == set(was)
    for domain in shown:
        if domain == "field-overlays":
            assert [field["createdAt"] for field in shown[domain]["userFields"]] == [None] * len(was[domain]["userFields"])
            continue
        assert tracker._hashes(shared.id)[domain] == golden.tracker()._hashes(golden.tracker().head_id)[domain], domain


def test_the_header_reads_as_written(golden: CellucidSession) -> None:
    session = EXPECTED["session"]
    assert golden.created_at == session["created_at"]
    assert (golden.app.version, golden.app.catalog) == (session["app"]["version"], session["app"]["catalog"])
    assert [stamp.type for stamp in golden.actions] == session["actions"]
    assert all(re.fullmatch(r"[0-9a-f]{64}", stamp.params) for stamp in golden.actions)
    scope = session["scope"]
    assert golden.scope.trackers == tuple(scope["trackers"])
    assert [(entry.tracker, entry.branches) for entry in golden.scope.branches] == [
        (entry["tracker"], tuple(entry["branches"])) for entry in scope["branches"]
    ]
    for key in ("content", "results", "notes", "chats", "journal", "times", "source"):
        assert getattr(golden.scope, key) == scope[key]
    assert (list(golden.requires.obs), list(golden.requires.var)) == (session["requires"]["obs"], session["requires"]["var"])
    # Page chrome outside the history (the pressed selection tool), read-only and never applied.
    assert golden.chrome == session["chrome"]
    golden.chrome["pressed"].clear()
    assert golden.chrome == session["chrome"]
    assert golden.active_tracker == session["active_tracker"] == T1
    assert golden.sha256 == GOLDEN_SHA256
    assert len(golden._contents.values) == session["values"]
    assert golden.dataset_fingerprint["cellOrder"]["digest"] == EXPECTED["dataset"]["cell_order_digest"]
    assert golden.source == {"kind": "sample", "dataset": EXPECTED["dataset"]["id"]}
    assert golden.cells is None and golden.genes is None


@pytest.mark.core
@pytest.mark.parametrize("index", [0, 1, 2], ids=["annotation", "replay", "removed"])
def test_every_tracker_reads_as_expected(golden: CellucidSession, index: int) -> None:
    expected = EXPECTED["trackers"][index]
    tracker = golden.trackers[index]
    assert (tracker.id, tracker.name, tracker.removed is not None, tracker.start.kind) == (
        expected["id"], expected["name"], expected["removed"], expected["start"]
    )
    assert (tracker.head_id, tracker.active_branch, tracker.main_branch) == (
        expected["head"], expected["active_branch"], expected["main_branch"]
    )
    paused = tracker.paused
    assert (None if paused is None else {"at": paused.at, "origins": [by.origin for by in paused.by], "draft": paused.draft}) == expected["paused"]
    assert [
        {"id": group.id, "name": group.name, "branch": group.branch, "first": group.first, "last": group.last,
         "steps": list(group.step_ids), "n_notes": len(group.notes)}
        for group in tracker.groups
    ] == expected["groups"]
    assert [
        {"step": step.id, "cause": step.squashed.cause, "count": step.squashed.count,
         "origins": [by.origin for by in step.squashed.by], "touches": list(step.touches),
         "carried": {"notes": list(step.squashed.carried.notes), "bookmarks": list(step.squashed.carried.bookmarks)}}
        for step in tracker.step_records if step.kind == "squashed"
    ] == expected["squashed"]
    assert [
        {"id": branch.id, "name": branch.name, "base": branch.base, "tip": branch.tip, "main": branch.main,
         "locked": branch.locked, "archived": branch.archived, "n_notes": len(branch.notes),
         "path": branch.steps["id"].tolist(), "own": branch.own_steps["id"].tolist()}
        for branch in tracker.branches
    ] == expected["branches"]
    table = tracker.steps
    assert list(table.columns) == list(STEP_COLUMNS) == expected["steps"]["columns"]
    rows = []
    for row in table.itertuples(index=False):
        values = list(row)
        at = values[STEP_COLUMNS.index("at")]
        values[STEP_COLUMNS.index("at")] = None if pd.isna(at) else int(at.value // 1_000_000)
        rows.append([None if value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value)) else value for value in values])
    assert rows == expected["steps"]["rows"]
    assert str(table["at"].dtype) == "datetime64[ns, UTC]"
    assert [
        {"id": turn.id, "label": turn.label, "origin": turn.by.origin, "branch": turn.branch,
         "steps": table.loc[table["turn"] == turn.id, "id"].tolist(), "summary": turn.notes[0].text if turn.notes else None}
        for turn in tracker.turns
    ] == expected["turns"]
    notes = tracker.notes_table()
    assert list(notes.columns) == list(NOTE_COLUMNS)
    assert [
        {"id": row.id, "on": {row.on: row.target}, "text": row.text, "origin": row.origin}
        for row in notes.itertuples(index=False)
    ] == expected["notes"]
    assert [{"id": entry.id, "op": entry.op, "label": entry.label, "purged": entry.purged} for entry in tracker.journal] == expected["journal"]
    assert tracker.stale() == expected["stale"]
    assert [entry["conversation"] for entry in golden.conversations(tracker)] == expected["conversations"]


def test_the_state_at_a_step_is_the_one_written(golden: CellucidSession) -> None:
    for state in EXPECTED["states"].values():
        tracker = golden.tracker(state["tracker"])
        assert tracker._hashes(state["step"]) == state["values"]
        values = tracker.step(state["step"]).state()
        assert set(values) == set(state["values"])
        assert all(_hash_value(golden._value(value_hash)) == value_hash for value_hash in state["values"].values())
    head = golden.tracker().step("s13").state()
    score = next(field for field in head["field-overlays"]["userFields"] if field["key"] == "score")
    assert score["values"]["cells"].tolist() == EXPECTED["values"]["score"]["cells"]
    assert score["values"]["values"].dtype == np.float32
    assert [None if math.isnan(value) else value for value in score["values"]["values"].tolist()] == pytest.approx(
        EXPECTED["values"]["score"]["values"], nan_ok=True
    )
    merged = next(field for field in head["field-overlays"]["userFields"] if field["key"] == "cell_type (merged)")
    assert merged["categories"] == EXPECTED["values"]["merged"]["categories"]
    assert merged["codes"]["values"].tolist() == EXPECTED["values"]["merged"]["codes"]
    assert head["highlights"]["pages"][0]["groups"][0]["cells"].tolist() == EXPECTED["values"]["doublets"]
    # Numbers JSON cannot hold come back as they were written.
    assert math.copysign(1.0, head["camera"]["views"][0]["orbit"]["theta"]) < 0
    assert math.isnan(head["settings"]["controls"]["gamma"]["value"])


def test_diffs_say_what_differs_domain_by_domain(golden: CellucidSession) -> None:
    for entry in EXPECTED["diffs"]:
        a = golden.tracker(entry["a"]["tracker"]).step(entry["a"]["step"])
        b = golden.tracker(entry["b"]["tracker"]).step(entry["b"]["step"])
        diff = a.diff(b)
        assert isinstance(diff, Diff)
        assert [{"domain": domain, "name": name, "line": line} for domain, name, line in diff.rows] == entry["rows"]
        assert list(diff.table.columns) == ["domain", "name", "line"]
        if entry["a"]["tracker"] == entry["b"]["tracker"]:
            assert diff.common is not None
        else:
            assert diff.common is None
    tracker = golden.tracker()
    assert tracker.step("s5").diff("s13").common == "s5"
    assert tracker.step("s13").diff("s9").common == "s2"


def test_methods_and_pseudocode_are_the_web_apps_words(golden: CellucidSession) -> None:
    for entry in EXPECTED["methods"]:
        assert golden.tracker(entry["tracker"]).methods(entry["branch"]) == entry["text"]
    for entry in EXPECTED["pseudocode"]:
        assert golden.tracker(entry["tracker"]).pseudocode(entry["branch"]) == entry["text"]
    # Main when no branch is named.
    assert golden.tracker().methods() == EXPECTED["methods"][0]["text"]


def test_records_are_found_by_name_or_id_and_refused_naming_the_others(golden: CellucidSession) -> None:
    # The removed tracker shares the active one's name: names are among the trackers not removed.
    assert golden.tracker("Annotation").id == T1
    assert golden.tracker(" annotation ").id == T1
    assert golden.tracker(T3).removed is not None
    with pytest.raises(ValueError) as refused:
        golden.tracker("Nothing")
    assert str(refused.value) == 'There is no tracker "Nothing"; the trackers are "Annotation" and "Annotation (replay)".'
    tracker = golden.tracker()
    assert tracker.branch("without SCORE").id == "b1"
    assert tracker.branch("b2").name == "Clusters"
    assert tracker.branch().id == tracker.active_branch
    with pytest.raises(ValueError) as refused:
        tracker.branch("Gone")
    assert str(refused.value) == 'There is no branch "Gone"; the branches are "Main", "Without score", "Clusters", "Scratch" and "Tidy fork".'
    with pytest.raises(ValueError) as refused:
        tracker.step("s11")
    # The reader names its own call, not the agents' history_get.
    assert str(refused.value) == 'There is no step "s11" in the tracker "Annotation"; Tracker.steps lists its steps.'


@pytest.mark.core
def test_steps_carry_every_field_of_the_model(golden: CellucidSession) -> None:
    tracker = golden.tracker()
    agent = tracker.step("s5")
    assert (agent.by.origin, agent.by.agent, agent.by.model, agent.by.turn) == ("agent", "Claude", "claude-sonnet-4-6", "claude-1")
    assert agent.bookmark == "Markers saved" and [note.text for note in agent.notes] == ["Top hits match the paper."]
    assert (agent.kind, agent.v, agent.group, agent.squashed) == ("action", 1, None, None)
    assert not hasattr(agent, "picture") and not hasattr(agent, "picture_png"), "no step holds a picture"
    root = tracker.step("s0")
    assert (root.kind, root.type, root.v) == ("root", None, None)
    # The notebook's whole turn, grouped.
    qc = tracker.group("qc")
    assert (qc.id, qc.name, qc.step_ids, qc.branch) == ("g0", "QC", ("s2", "s3"), "b0")
    assert tracker.group("g0") is qc and [note.text for note in qc.notes] == ["Doublets marked and the score pushed from the notebook."]
    assert tracker.step("s2").group == "g0" and qc.steps["id"].tolist() == ["s2", "s3"]
    with pytest.raises(ValueError) as refused:
        tracker.group("Clustering")
    assert str(refused.value) == 'There is no group "Clustering"; the groups are "QC".'
    # What changed while recording was paused: one squashed step with no action of its own.
    paused = tracker.step("s14")
    assert (paused.kind, paused.type, paused.params, paused.v, paused.code) == ("squashed", None, None, None, None)
    assert (paused.squashed.cause, paused.squashed.count, [by.origin for by in paused.squashed.by]) == ("paused", None, ["user", "host"])
    assert paused.label == "While paused: Colors, Settings" and paused.touches == ("colors", "settings")
    assert paused.state()["settings"]["controls"]["point-size"]["value"] == 5, "a squashed step restores exactly"
    paused_tracker = golden.tracker(T2).paused
    assert paused_tracker is not None and [by.origin for by in paused_tracker.by] == ["user"] and paused_tracker.draft is None
    # A squash in mid-branch keeps the last replaced step's id; what the replaced steps carried is kept.
    tidied = tracker.step("s16")
    assert (tidied.kind, tidied.squashed.cause, tidied.squashed.count) == ("squashed", "steps", 2)
    assert list(tidied.squashed.carried.bookmarks) == ["Before tidying"] and len(tidied.squashed.carried.notes) == 1
    assert tracker.step("s17").parent == "s16" and tracker.branch("Tidy fork").base == "s16"
    with pytest.raises(ValueError):
        tracker.step("s15")
    analysis = tracker.step("s4")
    assert analysis.stale == {"cause": "s13", "domains": ["field-overlays"]}
    assert set(analysis.reads) == {"cells", "mapping", "expression", "field-overlays", "views", "filters", "visibility", "highlights"}
    camera = tracker.step("s10")
    assert (camera.merged, camera.until - camera.at, camera.ms, camera.coalesce_key) == (3, 2000, None, "camera:View 1")
    flagged = tracker.step("s12")
    assert flagged.status == "unapplicable" and flagged.before == flagged.after
    assert flagged.reason == 'There is no highlight "Gate" on page "QC"; its highlights are "Doublets".'
    pushed = tracker.step("s3")
    assert pushed.params["cells"].dtype == np.uint32 and pushed.params["values"].dtype == np.float32
    replayed = golden.tracker(T2).step("s1")
    assert (replayed.replay_of.tracker, replayed.replay_of.step) == (T1, "s1")
    start = golden.tracker(T2).start
    assert start.from_.branch.name == "Main" and [entry.outcome for entry in start.report] == ["applied", "applied", "unmatched"]
    assert start.mapping.genes[0].from_ == "Cd3e"
    assert [entry.op for entry in tracker.journal] == [
        "note.add", "bookmark.set", "group.create", "branch.lock", "branch.archive", "step.delete", "squash",
    ]
    grouped = tracker.journal[2].after
    # A group edit records which steps the group held (memberships), never the steps.
    assert [group.name for group in grouped.groups] == ["QC"] and grouped.steps == ()
    assert [(item.group, item.from_, item.to) for item in grouped.memberships] == [("g0", "s2", "s3")]
    deleted = tracker.journal[5]
    assert deleted.op == "step.delete" and [step.id for step in deleted.before.steps] == ["s6", "s11", "s12"]
    # A step in a journal entry is a record of what was removed, not a step of the tree.
    assert type(deleted.before.steps[1]) is cellucid._history_model.Step


# =============================================================================
# Applying the golden fixture: the shared plans
# =============================================================================


@pytest.mark.core
def test_the_apply_sequence_gives_the_shared_plans(golden: CellucidSession) -> None:
    adata = golden_adata()
    assert cell_order_digest(adata.obsm["X_umap"]) == EXPECTED["dataset"]["cell_order_digest"]
    first = None
    for entry in EXPECTED["apply"]["sequence"]:
        plan = golden.apply(adata, tracker=entry["tracker"], branch=entry["branch"], step=entry["step"])
        assert isinstance(plan, ApplyPlan)
        assert (plan.tracker, plan.branch, plan.step) == (entry["expect"]["tracker"], entry["expect"]["branch"], entry["expect"]["step"])
        assert plan.session_sha256 == GOLDEN_SHA256
        assert {"obs": list(plan.requires.obs), "var": list(plan.requires.var)} == entry["requires"]
        assert _rows(plan) == entry["rows"], entry["name"]
        assert list(plan.table.columns) == ["kind", "domain", "name", "column", "from", "cells", "status", "note"]
        before = {name: adata.var[name].tolist() for name in adata.var.columns}
        report = plan.apply()
        assert isinstance(report, ApplyReport) and report.adata is adata
        # The var columns the apply wrote: new ones, and replaced ones (their values changed).
        written = [
            {"name": name, "values": adata.var[name].tolist()}
            for name in adata.var.columns if before.get(name) != adata.var[name].tolist()
        ]
        assert written == entry["var_columns"], entry["name"]
        assert all(adata.var[column["name"]].dtype == bool for column in written)
        assert report.applied is entry["applied"], entry["name"]
        assert (report.provenance is None) is not entry["applied"]
        if first is None:
            first = report
    assert first is not None and first.provenance is not None
    provenance = cellucid._history_model.as_json("Provenance", first.provenance)
    assert {key: value for key, value in provenance.items() if key not in ("applied_at", "cellucid_version", "package")} == EXPECTED["apply"]["provenance"]
    assert provenance["package"] == "python" and provenance["cellucid_version"] == package_version()
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z", provenance["applied_at"])
    # What the object holds now.
    assert adata.obs["cell_type"].cat.categories.tolist() == ["B cell", "T cell", "NK"]
    assert "cluster" in adata.obs and "leiden" not in adata.obs
    assert adata.obs["cellucid_highlight__QC__Doublets"].sum() == len(EXPECTED["values"]["doublets"])
    assert adata.obs["cellucid_highlight__Types__cell_type_T_cell"].sum() == len(EXPECTED["values"]["t_cells"])
    assert adata.obs["score"].dtype == np.float32
    assert adata.obs["cell_type (merged)"].cat.categories.tolist() == EXPECTED["values"]["merged"]["categories"]
    gene_set = adata.uns["cellucid"]["gene_sets"]["T cell markers"]
    assert (gene_set["genes"], gene_set["origin"]) == (["CD3E", "NKG7"], "session")
    marked = adata.var["cellucid_gene_set__T_cell_markers"]
    assert marked.dtype == bool
    assert adata.var_names[marked.to_numpy()].tolist() == ["CD3E", "NKG7"]
    records = history_applied(adata)
    assert [(record["step"], record["branch"]) for record in records] == [("s13", "b0"), ("s8", "b1")]
    assert json.loads(adata.uns["cellucid"]["history_applied"]) == records
    assert cellucid.history_applied(adata) == records


def test_applying_again_is_a_no_op_and_dry_run_false_applies_at_once(golden: CellucidSession) -> None:
    adata = golden_adata()
    report = apply_cellucid_session_to_anndata(adata, GOLDEN, dry_run=False)
    assert isinstance(report, ApplyReport) and report.applied
    again = apply_cellucid_session_to_anndata(adata, golden, dry_run=False)
    assert not again.applied and again.provenance is None
    assert {row.status for row in again.plan.rows} == {"unchanged", "not-applied"}
    assert len(history_applied(adata)) == 1


def test_keep_session_keeps_the_file_to_rebuild_the_history(golden: CellucidSession, tmp_path: Path) -> None:
    adata = golden_adata()
    report = golden.apply(adata, dry_run=False, keep_session=True)
    assert report.provenance is not None and report.provenance.session is not None
    rebuilt = tmp_path / "rebuilt.cellucid-session"
    rebuilt.write_bytes(base64.b64decode(report.provenance.session))
    assert rebuilt.read_bytes() == GOLDEN.read_bytes()
    assert read_session(rebuilt).tracker().steps["id"].tolist() == golden.tracker().steps["id"].tolist()
    # The record is JSON text, which .h5ad and Zarr hold as it is.
    path = tmp_path / "applied.h5ad"
    adata.write_h5ad(path)
    assert history_applied(ad.read_h5ad(path)) == history_applied(adata)


def test_keep_session_reads_the_file_when_the_plan_is_made(tmp_path: Path) -> None:
    """A plan made with keep_session holds the file's bytes: the file may go
    before the plan is applied (a viewer's session file is removed once
    planned), and a file changed since it was read is refused."""
    copy_path = tmp_path / "copy.cellucid-session"
    copy_path.write_bytes(GOLDEN.read_bytes())
    session = read_session(copy_path)
    plan = session.apply(golden_adata(), keep_session=True)
    copy_path.unlink()
    report = plan.apply()
    assert report.provenance is not None and report.provenance.session is not None
    assert base64.b64decode(report.provenance.session) == GOLDEN.read_bytes()
    changed = tmp_path / "changed.cellucid-session"
    changed.write_bytes(GOLDEN.read_bytes())
    later = read_session(changed)
    changed.write_bytes(GOLDEN.read_bytes() + b"\x00")
    with pytest.raises(ValueError, match="changed.cellucid-session changed since it was read"):
        later.apply(golden_adata(), keep_session=True)
    # Without keep_session the file is not needed at all.
    assert later.apply(golden_adata()).apply().applied


def test_a_column_the_person_made_is_never_replaced(golden: CellucidSession) -> None:
    adata = golden_adata()
    adata.obs["score"] = np.arange(adata.n_obs, dtype=np.float64)
    plan = golden.apply(adata)
    row = next(row for row in plan.rows if row.name == "score")
    assert (row.status, row.note) == (
        "not-applied", 'The object already has an obs column "score" that Cellucid did not write; it is never replaced.'
    )
    report = plan.apply()
    assert report.applied and adata.obs["score"].tolist() == list(range(adata.n_obs))
    assert report.provenance is not None and "score" not in report.provenance.columns
    assert any(item.kind == "user-field" and item.name == "score" for item in report.provenance.skipped)


def test_an_earlier_applys_column_is_replaced_only_with_overwrite(tmp_path: Path) -> None:
    n = 6
    adata = _small(n)
    root = files.state()
    one = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [0, 1])))
    two = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [2, 3, 4])))
    path = files.session_file(
        tmp_path / "gate.cellucid-session", fingerprint_record=_fingerprint(n),
        trackers=[{"id": "t1", "name": "Tracker 1", "root": root, "steps": [{"set": {"highlights": one}}, {"set": {"highlights": two}}]}],
    )
    session = read_session(path)
    session.apply(adata, step="s1", dry_run=False)
    plan = session.apply(adata, step="s2")
    row = next(row for row in plan.rows if row.kind == "highlight")
    assert (row.status, row.note) == ("not-applied", 'An earlier apply (step s1) wrote "cellucid_highlight__QC__Gate" with other values; overwrite replaces it.')
    assert not plan.writes_anything
    replaced = session.apply(adata, step="s2", overwrite=True, dry_run=False)
    assert [row.status for row in replaced.plan.rows if row.kind == "highlight"] == ["replace"]
    assert adata.obs["cellucid_highlight__QC__Gate"].tolist() == [False, False, True, True, True, False]
    assert [record["step"] for record in history_applied(adata)] == ["s1", "s2"]


def test_a_plan_refuses_an_object_changed_since(golden: CellucidSession) -> None:
    adata = golden_adata()
    plan = golden.apply(adata)
    adata.obs["score"] = 1.0
    with pytest.raises(ValueError, match="The object changed since this plan was made"):
        plan.apply()
    assert "cellucid" not in adata.uns


def test_the_cells_must_be_the_sessions_in_its_order(golden: CellucidSession) -> None:
    adata = golden_adata()
    order = np.r_[1, 0, np.arange(2, adata.n_obs)]
    shuffled = adata[order].copy()
    with pytest.raises(ValueError, match=r"This object's cells are not the cells the session was saved on, in the same order"):
        golden.apply(shuffled)
    # Other coordinates (an embedding computed again) are other cells too: the digest is over the stored values.
    moved = golden_adata()
    moved.obsm["X_umap"] = moved.obsm["X_umap"] + 1
    with pytest.raises(ValueError, match="cell-order digest of obsm\\['X_umap'\\] \\(columns 0, 1\\)"):
        golden.apply(moved)


def test_the_dataset_must_be_the_sessions(golden: CellucidSession) -> None:
    adata = golden_adata()
    adata.uns["cellucid"] = {"id": "another"}
    with pytest.raises(ValueError, match="saved on the dataset 'history_golden', and this object is 'another'"):
        golden.apply(adata)
    with pytest.raises(ValueError, match="saved on 600 cells, and this object has 599"):
        golden.apply(golden_adata()[:599].copy())


def test_only_what_the_chosen_step_uses_is_required(golden: CellucidSession) -> None:
    adata = golden_adata()
    adata.obs = adata.obs.rename(columns={"cell_type": "kind"})
    # The bookmark colors by cell_type; the head colors by a gene.
    with pytest.raises(ValueError, match="The state at step s5 uses the obs column 'cell_type'"):
        golden.apply(adata, step="s5")
    head = golden.apply(adata)
    assert [row.status for row in head.rows if row.kind == "categories-renamed"] == ["not-applied"]
    assert next(row.note for row in head.rows if row.kind == "categories-renamed") == 'The object has no obs column "cell_type".'


def test_a_step_off_the_branch_named_is_refused(golden: CellucidSession) -> None:
    with pytest.raises(ValueError, match=r'Step s13 is not on the branch "Clusters"'):
        golden.apply(golden_adata(), branch="Clusters", step="s13")
    plan = golden.apply(golden_adata(), branch="Clusters", step="s2")
    assert (plan.branch, plan.step) == ("b2", "s2")


# =============================================================================
# Apply rules, one by one (synthesized files)
# =============================================================================


def _small(n: int = 6) -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(["B", "T", "T", "NK", "B", "T"][:n], categories=["B", "T", "NK"]),
            "leiden": pd.Categorical([str(i % 3) for i in range(n)]),
            "label": np.array([f"x{i % 2}" for i in range(n)], dtype=object),
            "count": np.arange(n, dtype=np.int64),
        },
        index=pd.Index([f"cell_{i}" for i in range(n)]),
    )
    adata = ad.AnnData(X=np.ones((n, 2), dtype=np.float32), obs=obs, var=pd.DataFrame(index=pd.Index(["G1", "G2"])))
    adata.obsm["X_umap"] = np.column_stack([np.arange(n), np.arange(n) * 2]).astype(np.float32)
    return adata


def _fingerprint(n: int = 6, kept: np.ndarray | None = None) -> dict[str, Any]:
    """The fingerprint of ``_small(n)`` (its cell-order digest over ``kept``, or every cell)."""
    values = _small(n).obsm["X_umap"]
    digest = cell_order_digest(values if kept is None else values[kept], kept)
    return files.fingerprint(cell_count=n if kept is None else int(kept.size), var_count=2, dataset_id="exact", digest=digest)


def _write(tmp_path: Path, *steps: dict[str, Any], root: dict[str, Any] | None = None, n: int = 6, **options: Any) -> Path:
    return files.session_file(
        tmp_path / "small.cellucid-session",
        fingerprint_record=_fingerprint(n),
        trackers=[{"id": "t1", "name": "Tracker 1", "root": root or files.state(), "steps": list(steps)}],
        **options,
    )


def test_labels_rename_columns_and_categories_in_place(tmp_path: Path) -> None:
    overlays = files.overlays(
        renamed_fields=[{"source": "obs", "key": "leiden", "name": "cluster"}, {"source": "obs", "key": "absent", "name": "x"},
                        {"source": "var", "key": "G1", "name": "Gene one"}],
        renamed_categories=[
            {"source": "obs", "field": "cell_type", "category": "B", "name": "B cell"},
            {"source": "obs", "field": "cell_type", "category": "NK", "name": "T"},
            {"source": "obs", "field": "label", "category": "x1", "name": "y"},
            {"source": "obs", "field": "count", "category": 1, "name": "one"},
            {"source": "obs", "field": "leiden", "category": 1, "name": "One"},
        ],
        deleted_fields=[{"source": "obs", "key": "label", "purged": False}],
    )
    adata = _small()
    plan = read_session(_write(tmp_path, {"set": {"field-overlays": overlays}})).apply(adata)
    assert [(row.kind, row.name, row.column, list(row.from_), row.status, row.note) for row in plan.rows] == [
        ("column-renamed", "cluster", None, ["leiden"], "write", None),
        ("column-renamed", "x", None, ["absent"], "not-applied", 'The object has no obs column "absent".'),
        ("column-renamed", "Gene one", None, ["G1"], "not-applied", "A gene's display name is recorded, never applied to the object."),
        ("categories-renamed", "B cell", "cell_type", ["B"], "write", None),
        ("categories-merged", "T", "cell_type", ["NK"], "write", None),
        ("categories-renamed", "y", "label", ["x1"], "write", None),
        ("categories-renamed", "one", "count", ["1"], "not-applied", 'The obs column "count" holds numbers, not categories.'),
        # leiden holds "1" (text), not the number 1: labels told apart as JavaScript tells them.
        ("categories-renamed", "One", "cluster", ["1"], "not-applied", 'The obs column "cluster" has no category "1".'),
        ("field-deleted", "label", None, [], "not-applied", "A deletion is recorded, never applied to the object: the column stays."),
    ]
    assert [row.cells for row in plan.rows if row.kind.startswith("categories-") and row.status == "write"] == [2, 1, 3]
    plan.apply()
    assert list(adata.obs.columns) == ["cell_type", "cluster", "label", "count"]
    assert adata.obs["cell_type"].cat.categories.tolist() == ["B cell", "T"]
    assert adata.obs["cell_type"].tolist() == ["B cell", "T", "T", "T", "B cell", "T"]
    assert adata.obs["label"].tolist() == ["x0", "y", "x0", "y", "x0", "y"]


def test_labels_of_boolean_and_identifier_columns_are_the_viewers(tmp_path: Path) -> None:
    """A boolean column's categories are false and true, and an integer column
    holding identifiers past 2^53 has its exact decimal texts as categories,
    as the viewer reads them: a rename made in the viewer applies."""
    first = 2**60
    overlays = files.overlays(renamed_categories=[
        {"source": "obs", "field": "flag", "category": True, "name": "kept"},
        {"source": "obs", "field": "donor_id", "category": str(first + 1), "name": "donor B"},
        {"source": "obs", "field": "count", "category": 1, "name": "one"},
    ])
    adata = _small()
    adata.obs["flag"] = np.array([True, False, True, False, False, True])
    adata.obs["donor_id"] = np.array([first + i for i in range(6)], dtype=np.int64)
    plan = read_session(_write(tmp_path, {"set": {"field-overlays": overlays}})).apply(adata)
    assert [(row.kind, row.name, row.column, list(row.from_), row.status) for row in plan.rows] == [
        ("categories-renamed", "kept", "flag", ["true"], "write"),
        ("categories-renamed", "donor B", "donor_id", [str(first + 1)], "write"),
        # Integers all within 2^53 stay numbers.
        ("categories-renamed", "one", "count", ["1"], "not-applied"),
    ]
    assert [row.cells for row in plan.rows if row.status == "write"] == [3, 1]
    plan.apply()
    assert adata.obs["flag"].tolist() == ["kept", "false", "kept", "false", "false", "kept"]
    assert adata.obs["donor_id"].tolist() == [str(first), "donor B", *[str(first + i) for i in range(2, 6)]]


def test_swapped_column_names_rename_at_once(tmp_path: Path) -> None:
    overlays = files.overlays(renamed_fields=[
        {"source": "obs", "key": "leiden", "name": "label"}, {"source": "obs", "key": "label", "name": "leiden"},
        {"source": "obs", "key": "count", "name": "cell_type"},
    ])
    adata = _small()
    plan = read_session(_write(tmp_path, {"set": {"field-overlays": overlays}})).apply(adata, dry_run=False).plan
    assert [row.status for row in plan.rows] == ["write", "write", "not-applied"]
    assert adata.obs["label"].tolist() == _small().obs["leiden"].tolist()


def test_user_fields_apply_by_kind(tmp_path: Path) -> None:
    fields = files.overlays(user_fields=[
        files.user_category("user_cat_1", "annotation", ["a", "b"], [0, 1, 255, 0, 1, 0]),
        files.user_values("user_obs_value_1", "pushed", [0.5, math.nan, 2.0], cells_of=[0, 2, 4]),
        files.user_copy("user_obs_value_2", "copy of count", "obs", "count"),
        files.user_category("user_cat_2", "deleted one", ["c"], [0] * 6, deleted=True),
    ])
    adata = _small()
    plan = read_session(_write(tmp_path, {"set": {"field-overlays": fields}})).apply(adata)
    assert [(row.kind, row.name, row.cells, row.status, row.note) for row in plan.rows] == [
        ("user-field", "annotation", 5, "write", None),
        ("user-field", "pushed", 2, "write", None),
        ("user-field", "copy of count", 6, "write", None),
        ("user-field", "deleted one", None, "not-applied", "Deleted in the session."),
    ]
    plan.apply()
    assert adata.obs["annotation"].astype(object).where(adata.obs["annotation"].notna(), None).tolist() == ["a", "b", None, "a", "b", "a"]
    assert adata.obs["pushed"].dtype == np.float32
    assert np.isnan(adata.obs["pushed"].to_numpy()[[1, 2, 3, 5]]).all()
    assert adata.obs["copy of count"].tolist() == [0, 1, 2, 3, 4, 5]
    assert "deleted one" not in adata.obs


def test_a_user_category_fields_palette_is_read_and_checked(tmp_path: Path) -> None:
    # The colors a field made from highlight pages opens with (its pages'), one
    # #rrggbb a category, or null.
    made = {**files.user_category("user_cat_1", "annotation", ["a", "b"], [0, 1, 255, 0, 1, 0]),
            "palette": ["#dc2626", "#16a34a"]}
    plan = read_session(_write(tmp_path, {"set": {"field-overlays": files.overlays(user_fields=[made])}})).apply(_small())
    assert [row.status for row in plan.rows] == ["write"]
    for palette in (["#dc2626"], ["#DC2626", "#16a34a"], ["red", "#16a34a"], "#dc2626"):
        bad = {**made, "palette": palette}
        with pytest.raises(ValueError, match=r"palette must be null or one #rrggbb color a category"):
            read_session(_write(tmp_path, {"set": {"field-overlays": files.overlays(user_fields=[bad])}})).apply(_small())


def test_a_field_of_the_intersections_of_41_pages_applies(tmp_path: Path) -> None:
    # Each intersection combines any number of pages; the session keys
    # its label by its pages' positions ("0+1+...+40"), a record of labels;
    # the column is its codes.
    pages = [{"pageId": f"page_{p}", "label": f"P{p}"} for p in range(41)]
    every = "+".join(str(p) for p in range(41))
    categories = [f"P{p}" for p in range(41)] + ["First two", "Every page"]
    for labels in ({"0+1": "First two", every: "Every page"},):
        made = {**files.user_category("user_cat_1", "pages", categories, [42, 41, 0, 1, 2, 40]),
                "sourcePages": pages, "overlapStrategy": "intersections", "intersectionLabels": labels}
        adata = _small()
        plan = read_session(_write(tmp_path, {"set": {"field-overlays": files.overlays(user_fields=[made])}})).apply(adata)
        assert [row.status for row in plan.rows] == ["write"]
        plan.apply()
        assert adata.obs["pages"].tolist() == ["Every page", "First two", "P0", "P1", "P2", "P40"]


def test_a_written_columns_colors_are_the_ones_it_is_shown_in(tmp_path: Path) -> None:
    # A field's palette goes with its column as Scanpy keeps colors,
    # uns['<key>_colors'], one a category in the column's order -- every
    # category, a used or not. A field without a palette opens in
    # Cellucid's own category colors, and is written with them; a
    # column of no categories (a highlight) takes no colors: those under its
    # name described another column.
    made = {**files.user_category("user_cat_1", "annotation", ["a", "b", "unused"], [0, 1, 255, 0, 1, 0]),
            "palette": ["#dc2626", "#16a34a", "#2563eb"]}
    plain = files.user_category("user_cat_2", "plain", ["x"], [0] * 6)
    groups = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [0, 1])))
    path = _write(tmp_path, {"set": {"field-overlays": files.overlays(user_fields=[made, plain]), "highlights": groups}})
    adata = _small()
    adata.uns["plain_colors"] = ["#000000", "#111111", "#222222"]
    adata.uns["cellucid_highlight__QC__Gate_colors"] = np.array(["#000000"], dtype=object)
    adata.uns["cell_type_colors"] = ["#aaaaaa", "#bbbbbb", "#cccccc"]
    plan = read_session(path).apply(adata)
    assert "annotation_colors" not in adata.uns and "plain_colors" in adata.uns  # a dry run writes nothing
    plan.apply()
    assert adata.obs["annotation"].cat.categories.tolist() == ["a", "b", "unused"]
    assert adata.uns["annotation_colors"] == ["#dc2626", "#16a34a", "#2563eb"]
    assert adata.uns["plain_colors"] == ["#ff9900"]
    assert "cellucid_highlight__QC__Gate_colors" not in adata.uns
    assert adata.uns["cell_type_colors"] == ["#aaaaaa", "#bbbbbb", "#cccccc"]  # a column left as it is keeps its colors
    # Applied again: unchanged, colors as they are.
    adata.uns["annotation_colors"] = ["#dc2626", "#16a34a", "#000000"]
    assert not read_session(path).apply(adata, dry_run=False).applied
    assert adata.uns["annotation_colors"] == ["#dc2626", "#16a34a", "#000000"]
    # Replaced (overwrite): the new palette replaces the colors.
    other = {**made, "codes": {"cells": None, "values": files.typed("Uint8Array", [1] * 6)}, "palette": ["#111111", "#222222", "#333333"]}
    again = _write(tmp_path, {"set": {"field-overlays": files.overlays(user_fields=[other, plain]), "highlights": groups}})
    report = read_session(again).apply(adata, dry_run=False, overwrite=True)
    assert [row.status for row in report.plan.rows if row.name == "annotation"] == ["replace"]
    assert adata.uns["annotation_colors"] == ["#111111", "#222222", "#333333"]


def test_colors_follow_a_renamed_column_and_merged_categories(tmp_path: Path) -> None:
    overlays = files.overlays(
        renamed_fields=[{"source": "obs", "key": "leiden", "name": "cluster"}, {"source": "obs", "key": "count", "name": "number"}],
        renamed_categories=[
            # B merges into T, in B's place: the category keeps B's color.
            {"source": "obs", "field": "cell_type", "category": "B", "name": "T"},
            {"source": "obs", "field": "label", "category": "x1", "name": "y"},
        ],
    )
    adata = _small()
    adata.uns["cell_type_colors"] = np.array(["#aa0000", "#00aa00", "#0000aa"], dtype=object)
    adata.uns["leiden_colors"] = ["#111111", "#222222", "#333333"]
    adata.uns["cluster_colors"] = ["#999999"]  # stale: no column of that name
    adata.uns["count_colors"] = np.zeros((2, 4))  # not Scanpy's form: not carried
    adata.uns["label_colors"] = ["#444444", "#555555"]
    read_session(_write(tmp_path, {"set": {"field-overlays": overlays}})).apply(adata, dry_run=False)
    assert adata.obs["cell_type"].cat.categories.tolist() == ["T", "NK"]
    assert adata.uns["cell_type_colors"] == ["#aa0000", "#0000aa"]
    assert adata.uns["cluster_colors"] == ["#111111", "#222222", "#333333"]
    assert "leiden_colors" not in adata.uns
    assert "count_colors" not in adata.uns and "number_colors" not in adata.uns
    # A text column relabeled: Scanpy colors it by its sorted labels, which the rename reorders.
    assert "label_colors" not in adata.uns


def test_merged_categories_drop_colors_that_did_not_cover_them(tmp_path: Path) -> None:
    overlays = files.overlays(renamed_categories=[
        {"source": "obs", "field": "cell_type", "category": "NK", "name": "T"},
        {"source": "obs", "field": "leiden", "category": "0", "name": "zero"},
    ])
    adata = _small()
    adata.uns["cell_type_colors"] = ["#aa0000", "#00aa00"]  # fewer than its 3 categories: Scanpy draws its own
    adata.uns["leiden_colors"] = ["#111111"]  # a rename alone: left as it is
    read_session(_write(tmp_path, {"set": {"field-overlays": overlays}})).apply(adata, dry_run=False)
    assert adata.obs["cell_type"].cat.categories.tolist() == ["B", "T"]
    assert "cell_type_colors" not in adata.uns
    assert adata.uns["leiden_colors"] == ["#111111"]


def _recolor(field_id: str, *pairs: tuple[Any, str]) -> dict[str, Any]:
    """One field's entry of a view's colors: its categories recolored."""
    return {"field": {"source": "obs", "user": field_id}, "colors": [{"category": c, "color": h} for c, h in pairs]}


def _views(focus: int, count: int) -> dict[str, Any]:
    view = {"embedding": {"key": "X_umap", "columns": [0, 1]}, "followsDefault": True, "centroids": {"points": False, "labels": False}}
    return {"layout": "grid", "focus": focus, "liveHidden": False, "camerasLinked": True, "views": [view] * count}


def test_a_fields_colors_are_the_ones_the_view_in_focus_shows(tmp_path: Path) -> None:
    """The colors a person or an agent gave a field's categories reach
    uns['<key>_colors']: in the view in focus at the step, each
    category recolored there, the rest in the colors it opens with -- its
    palette, else Cellucid's own colors by position. Another view's colors
    stay how the data was shown."""
    stage = files.user_category("user_cat_1", "stage", ["early", "mid", "late"], [0, 1, 2, 0, 1, 2])
    fate = {**files.user_category("user_cat_2", "fate", ["a", "b"], [0, 1, 0, 1, 0, 1]), "palette": ["#dc2626", "#16a34a"]}
    live = [_recolor("user_cat_1", ("early", "#111111"))]
    kept = [
        {"field": {"source": "obs", "key": "cell_type"}, "colors": [{"category": "B", "color": "#000000"}]},
        _recolor("user_cat_1", ("late", "#008080"), ("early", "#808080"), ("absent", "#ffffff")),
        _recolor("user_cat_2", ("b", "#ffd700")),
    ]
    made = {"set": {"field-overlays": files.overlays(user_fields=[stage, fate])}}
    root = files.state(colors={"views": [[], []]}, views=_views(0, 2))

    def applied(focus: int) -> tuple[ad.AnnData, ApplyReport]:
        path = _write(tmp_path, made, {"set": {"colors": {"views": [live, kept]}, "views": _views(focus, 2)}}, root=root)
        adata = _small()
        return adata, read_session(path).apply(adata, dry_run=False)

    adata, report = applied(1)
    assert adata.uns["stage_colors"] == ["#808080", "#0077c8", "#008080"]
    assert adata.uns["fate_colors"] == ["#dc2626", "#ffd700"]
    # The dataset column recolored in the view in focus.
    assert adata.uns["cell_type_colors"] == ["#000000", "#0077c8", "#ad50e0"]
    column_row, colors_row = [row for row in report.plan.rows if row.domain == "colors"]
    assert (column_row.kind, column_row.name, column_row.status) == ("column-colors", "cell_type", "write")
    assert (colors_row.status, colors_row.note) == ("not-applied", COLORS_NOTE)
    assert next(row.note for row in report.plan.rows if row.domain == "views") == PRESENTATION_NOTE
    adata, _ = applied(0)
    assert adata.uns["stage_colors"] == ["#111111", "#0077c8", "#ad50e0"]
    assert adata.uns["fate_colors"] == ["#dc2626", "#16a34a"]
    assert "cell_type_colors" not in adata.uns  # recolored in another view only
    # The colors at the step before (none recolored): the colors the fields open with.
    path = _write(tmp_path, made, {"set": {"colors": {"views": [live, kept]}, "views": _views(1, 2)}}, root=root)
    adata = _small()
    read_session(path).apply(adata, step="s1", dry_run=False)
    assert adata.uns["stage_colors"] == ["#ff9900", "#0077c8", "#ad50e0"]
    # At the head: the columns are unchanged (their colors as they are); the
    # dataset column recolored there takes its colors, once.
    report = read_session(path).apply(adata, dry_run=False)
    assert [(row.kind, row.status) for row in report.plan.rows if row.domain == "colors"] == [
        ("column-colors", "write"), ("presentation", "not-applied")]
    assert adata.uns["stage_colors"] == ["#ff9900", "#0077c8", "#ad50e0"]
    report = read_session(path).apply(adata, dry_run=False)
    assert not report.applied
    assert [row.status for row in report.plan.rows if row.domain == "colors"] == ["unchanged", "not-applied"]
    assert next(row.note for row in report.plan.rows if row.domain == "colors" and row.kind == "presentation") == COLORS_NOTE


def _dataset_recolor(key: str, *pairs: tuple[Any, str]) -> dict[str, Any]:
    """A dataset column's entry of a view's colors: its categories recolored."""
    return {"field": {"source": "obs", "key": key}, "colors": [{"category": c, "color": h} for c, h in pairs]}


def test_a_recolored_dataset_column_is_written_with_the_colors_it_is_shown_in(tmp_path: Path) -> None:
    """A column that came with the dataset, recolored in the view in focus,
    takes the colors it is shown in: the categories recolored, the
    rest in the colors it opens with -- Scanpy's colors for it, as drawn
    (#rgb widened, alpha dropped), else Cellucid's own by position. Applied
    again: unchanged. A column it cannot be written to says why (a column
    the object lacks is refused before, as every field the state uses)."""
    recolors = [
        _dataset_recolor("cell_type", ("T", "#123456")),
        _dataset_recolor("leiden", ("1", "#00ff00")),
        _dataset_recolor("label", ("x0", "#ffffff")),
        _dataset_recolor("count", ("absent", "#ffffff")),
    ]
    path = _write(tmp_path, {"set": {"colors": {"views": [recolors]}, "views": _views(0, 1)}},
                  root=files.state(colors={"views": [[]]}, views=_views(0, 1)))
    adata = _small()
    adata.uns["cell_type_colors"] = np.array(["#AA0000", "#0a0", "#0000aaff"], dtype=object)
    adata.uns["leiden_colors"] = ["#111111"]  # fewer than its 3 categories: the viewer opens it in Cellucid's colors
    plan = read_session(path).apply(adata)
    rows = [(row.kind, row.name, row.status, row.note) for row in plan.rows if row.domain == "colors"]
    assert rows == [
        ("column-colors", "cell_type", "write", None),
        ("column-colors", "leiden", "write", None),
        ("column-colors", "label", "not-applied", 'The obs column "label" is not categorical; Scanpy keeps colors for a categorical column only.'),
        ("column-colors", "count", "not-applied", 'The obs column "count" is not categorical; Scanpy keeps colors for a categorical column only.'),
        ("presentation", "Colors", "not-applied", COLORS_NOTE),
    ]
    assert list(adata.uns["cell_type_colors"]) == ["#AA0000", "#0a0", "#0000aaff"]  # a dry run writes nothing
    report = plan.apply()
    assert report.applied and report.provenance is not None
    assert adata.uns["cell_type_colors"] == ["#aa0000", "#123456", "#0000aa"]
    assert adata.uns["leiden_colors"] == ["#ff9900", "#00ff00", "#ad50e0"]
    assert "label_colors" not in adata.uns
    assert [entry.name for entry in report.provenance.skipped if entry.kind == "column-colors"] == ["label", "count"]
    again = read_session(path).apply(adata, dry_run=False)
    assert not again.applied
    assert [row.status for row in again.plan.rows if row.kind == "column-colors"] == [
        "unchanged", "unchanged", "not-applied", "not-applied"]


def test_a_recolored_dataset_column_renamed_in_the_session_takes_its_colors_by_its_new_labels(tmp_path: Path) -> None:
    """The colors value names a dataset column by the key the dataset gives it
    and its categories by the labels shown at the step: the column renamed
    and a category relabeled by the same apply are written under the new
    names, the colors in the column's order; a label the column lacks is
    said, the others written."""
    overlays = files.overlays(
        renamed_fields=[{"source": "obs", "key": "cell_type", "name": "kind"}],
        renamed_categories=[{"source": "obs", "field": "cell_type", "category": "B", "name": "Bcell"}],
    )
    recolors = [_dataset_recolor("cell_type", ("Bcell", "#101010"), ("ghost", "#202020"))]
    path = _write(tmp_path, {"set": {"field-overlays": overlays, "colors": {"views": [recolors]}, "views": _views(0, 1)}},
                  root=files.state(colors={"views": [[]]}, views=_views(0, 1)))
    adata = _small()
    adata.uns["cell_type_colors"] = ["#aa0000", "#00aa00", "#0000aa"]
    report = read_session(path).apply(adata, dry_run=False)
    assert adata.obs["kind"].cat.categories.tolist() == ["Bcell", "T", "NK"]
    assert adata.uns["kind_colors"] == ["#101010", "#00aa00", "#0000aa"]
    assert "cell_type_colors" not in adata.uns
    (row,) = [row for row in report.plan.rows if row.kind == "column-colors"]
    assert (row.name, row.status, row.note) == (
        "kind", "write", 'Not written: the color of "ghost", which the obs column "kind" does not have.')
    # Only labels the column lacks: nothing to write.
    only = [_dataset_recolor("leiden", ("9", "#202020"), ("8", "#303030"))]
    path = _write(tmp_path, {"set": {"colors": {"views": [only]}, "views": _views(0, 1)}},
                  root=files.state(colors={"views": [[]]}, views=_views(0, 1)))
    plan = read_session(path).apply(_small())
    (row,) = [row for row in plan.rows if row.kind == "column-colors"]
    assert (row.status, row.note) == (
        "not-applied", 'Not written: the colors of "9", "8", which the obs column "leiden" does not have.')


def test_colors_the_session_holds_are_checked(tmp_path: Path) -> None:
    stage = files.user_category("user_cat_1", "stage", ["early"], [0] * 6)
    made = {"set": {"field-overlays": files.overlays(user_fields=[stage])}}
    for value, views, message in (
        ({"views": [[_recolor("user_cat_1", ("early", "#FFFFFF"))]]}, None, r"color must be one #rrggbb color"),
        ({"views": [[{"field": {"source": "obs", "user": ""}}]]}, None, r"field.user must be a non-empty text"),
        ({"views": [[]]}, _views(1, 2), r"colors value holds 1 view, and the view in focus is view 1"),
        ({"views": [[]]}, {**_views(0, 1), "focus": -1}, r"focus must be a non-negative whole number"),
    ):
        changes: dict[str, Any] = {"colors": value} if views is None else {"colors": value, "views": views}
        path = _write(tmp_path, made, {"set": changes}, root=files.state(colors={"views": [[]]}, views=_views(0, 1)))
        with pytest.raises(ValueError, match=message):
            read_session(path).apply(_small())


def test_a_session_saved_under_the_id_an_object_was_shown_under_applies_when_named(tmp_path: Path) -> None:
    """A viewer given dataset_id= serves the object under that id, and its
    session is saved under it: naming it applies the session to the
    object, which declares its own id; another id is refused."""
    record = files.fingerprint(cell_count=6, var_count=2, dataset_id="pancreas-notebook",
                               digest=_fingerprint()["cellOrder"]["digest"])
    gate = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [0, 1])))
    path = files.session_file(
        tmp_path / "shown.cellucid-session", fingerprint_record=record,
        trackers=[{"id": "t1", "name": "Tracker 1", "root": files.state(), "steps": [{"set": {"highlights": gate}}]}],
    )
    adata = _small()
    adata.uns["cellucid"] = {"id": "pancreas"}
    with pytest.raises(ValueError, match=(r"saved on the dataset 'pancreas-notebook', and this object is 'pancreas' .*"
                                          r"name that id: dataset_id='pancreas-notebook'")):
        apply_cellucid_session_to_anndata(adata, path)
    with pytest.raises(ValueError, match=r"saved on the dataset 'pancreas-notebook', and this object was shown as 'pancreas'"):
        apply_cellucid_session_to_anndata(adata, path, dataset_id="pancreas")
    with pytest.raises(ValueError, match="dataset_id"):
        apply_cellucid_session_to_anndata(adata, path, dataset_id="no/slash")
    report = read_session(path).apply(adata, dataset_id="pancreas-notebook", dry_run=False)
    assert report.applied and adata.obs["cellucid_highlight__QC__Gate"].tolist() == [True, True, False, False, False, False]
    # A plan made so is applied so (the plan is made again first).
    again = _small()
    again.uns["cellucid"] = {"id": "pancreas"}
    assert apply_cellucid_session_to_anndata(again, path, dataset_id="pancreas-notebook").apply().applied
    # A pull of the session file takes the same id.
    pulled = _small()
    pulled.uns["cellucid"] = {"id": "pancreas"}
    with pytest.raises(ValueError, match="this object is 'pancreas'"):
        cellucid.pull(pulled, path)
    assert list(cellucid.pull(pulled, path, dataset_id="pancreas-notebook").columns) == ["cellucid_highlight__highlight_1"]


def test_a_pull_of_a_session_file_writes_a_fields_palette_under_its_column(tmp_path: Path) -> None:
    made = {**files.user_category("user_cat_1", "annotation", ["a", "b"], [0, 1, 255, 0, 1, 0]), "palette": ["#dc2626", "#16a34a"]}
    path = _write(tmp_path, {"set": {"field-overlays": files.overlays(user_fields=[made])}})
    adata = _small()
    adata.uns["field__annotation_colors"] = ["#000000"]
    summary = cellucid.pull(adata, path, field_prefix="field__")
    assert list(summary.columns) == ["field__annotation"]
    assert adata.uns["field__annotation_colors"] == ["#dc2626", "#16a34a"]
    # Written beside a column of that name (on_conflict="suffix"): under the name written.
    adata = _small()
    adata.obs["annotation"] = 1.0
    adata.uns["annotation_colors"] = ["#000000"]
    summary = cellucid.pull(adata, path, on_conflict="suffix")
    (written,) = summary.columns
    assert written != "annotation"
    assert adata.uns[f"{written}_colors"] == ["#dc2626", "#16a34a"]
    assert adata.uns["annotation_colors"] == ["#000000"]


def test_highlight_groups_of_every_page_are_named_by_page_and_label(tmp_path: Path) -> None:
    groups = files.highlights(
        files.page("page_1", "QC / pass", files.group("highlight_1", "Gate", [0, 1]), files.group("highlight_2", "Gate!", [2])),
        files.page("page_2", "Types", files.group("highlight_3", "cell_type: T", [1, 2, 5], kind="category",
                                                   field={"source": "obs", "key": "cell_type"}, category="T")),
    )
    adata = _small()
    plan = read_session(_write(tmp_path, {"set": {"highlights": groups}})).apply(adata, dry_run=False).plan
    assert [(row.name, row.cells) for row in plan.rows if row.kind == "highlight"] == [
        ("cellucid_highlight__QC_pass__Gate", 2),
        ("cellucid_highlight__QC_pass__Gate_2", 1),
        ("cellucid_highlight__Types__cell_type_T", 3),
    ]
    assert adata.obs["cellucid_highlight__QC_pass__Gate"].dtype == bool
    assert {"obs": list(plan.requires.obs), "var": list(plan.requires.var)} == {"obs": ["cell_type"], "var": []}


def test_a_subsample_knows_only_its_cells(tmp_path: Path) -> None:
    n = 6
    record, kept = files.subsample_record(n, 4, seed=3)
    groups = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [int(kept[0])])))
    path = files.session_file(
        tmp_path / "subsample.cellucid-session",
        fingerprint_record=_fingerprint(n, kept=np.asarray(kept, dtype=np.int64)),
        trackers=[{"id": "t1", "name": "Tracker 1", "root": files.state(), "steps": [{"set": {"highlights": groups}}]}],
        cells_record=record,
    )
    adata = _small(n)
    session = read_session(path)
    assert session.cells is not None and session.cells["kept"].tolist() == kept.tolist()
    session.apply(adata, dry_run=False)
    column = adata.obs["cellucid_highlight__QC__Gate"]
    assert str(column.dtype) == "boolean"
    unknown = sorted(set(range(n)) - set(kept.tolist()))
    assert column.isna().to_numpy()[unknown].all() and not column.isna().to_numpy()[kept].any()
    with pytest.raises(ValueError, match="saved on 4 of 6 cells, and this object has 5"):
        session.apply(_small(5))


def test_gene_sets_are_written_once_and_never_over_the_persons(tmp_path: Path) -> None:
    sets = {"sets": [{"name": "Markers", "description": "d", "genes": ["G1"]}, {"name": "Mine", "description": "", "genes": ["G2"]}]}
    adata = _small()
    adata.uns["cellucid"] = {"gene_sets": {"Mine": {"genes": ["G1", "G2"], "description": "", "origin": "user"}}}
    plan = read_session(_write(tmp_path, {"set": {"gene-sets": sets}})).apply(adata, dry_run=False).plan
    assert [(row.name, row.status, row.note) for row in plan.rows] == [
        ("Markers", "write", None),
        ("Mine", "not-applied", 'The object already has a gene set "Mine" that no session wrote; it is never replaced.'),
    ]
    assert adata.uns["cellucid"]["gene_sets"]["Mine"]["genes"] == ["G1", "G2"]
    assert adata.uns["cellucid"]["gene_sets"]["Markers"]["origin"] == "session"


#: The golden head's gene set's var column.
T_CELL_MARKERS = "cellucid_gene_set__T_cell_markers"


def test_a_gene_sets_var_column_is_replaced_in_place_only_with_overwrite(golden: CellucidSession) -> None:
    # An earlier apply wrote the set with other genes, and its column after another var column.
    adata = golden_adata()
    adata.var["highly_variable"] = True
    golden.apply(adata, dry_run=False)
    assert list(adata.var.columns) == ["highly_variable", T_CELL_MARKERS]
    adata.uns["cellucid"]["gene_sets"]["T cell markers"]["genes"] = ["LYZ"]
    adata.var[T_CELL_MARKERS] = [False, False, False, True]
    kept = golden.apply(adata, dry_run=False)
    assert [row.status for row in kept.plan.rows if row.kind == "gene-set"] == ["not-applied"]
    assert not kept.applied and adata.var[T_CELL_MARKERS].tolist() == [False, False, False, True]
    replaced = golden.apply(adata, overwrite=True, dry_run=False)
    assert [row.status for row in replaced.plan.rows if row.kind == "gene-set"] == ["replace"]
    assert adata.var[T_CELL_MARKERS].tolist() == [True, False, True, False]
    assert adata.var[T_CELL_MARKERS].dtype == bool
    assert list(adata.var.columns) == ["highly_variable", T_CELL_MARKERS], "replaced in its place"
    # The provenance's columns stay obs columns.
    assert replaced.provenance is not None and T_CELL_MARKERS not in replaced.provenance.columns


def test_a_var_column_the_dataset_holds_beside_no_applied_set_is_its_own(golden: CellucidSession) -> None:
    adata = golden_adata()
    adata.var[T_CELL_MARKERS] = np.zeros(adata.n_vars, dtype=bool)
    plan = golden.apply(adata)
    assert [row.status for row in plan.rows if row.kind == "gene-set"] == ["write"]
    message = f'The gene set "T cell markers" has no var column: the var column "{T_CELL_MARKERS}" is the dataset\'s own.'
    with pytest.warns(UserWarning, match=re.escape(message)):
        plan.apply()
    assert adata.var[T_CELL_MARKERS].tolist() == [False] * adata.n_vars
    assert adata.uns["cellucid"]["gene_sets"]["T cell markers"]["genes"] == ["CD3E", "NKG7"]


def _gene_set_session(tmp_path: Path, genes: list[str], *, data_mapping: dict[str, Any] | None = None, var_count: int = 2) -> Path:
    mapping = files.mapping() if data_mapping is None else data_mapping
    record = files.fingerprint(cell_count=6, var_count=var_count, dataset_id="exact", digest=_fingerprint()["cellOrder"]["digest"],
                               data_mapping=mapping)
    sets = {"sets": [{"name": "Markers", "description": "", "genes": genes}]}
    return files.session_file(
        tmp_path / "sets.cellucid-session", fingerprint_record=record,
        trackers=[{"id": "t1", "name": "Tracker 1", "root": files.state(mapping=mapping), "steps": [{"set": {"gene-sets": sets}}]}],
    )


@pytest.mark.parametrize(("expression", "var_count"), [("raw", 3), (None, 0)])
def test_a_mapping_that_reads_no_var_genes_keeps_the_gene_sets_in_uns_only(
    tmp_path: Path, expression: str | None, var_count: int,
) -> None:
    adata = _small()
    adata.raw = ad.AnnData(X=np.ones((6, 3), dtype=np.float32), obs=adata.obs[[]], var=pd.DataFrame(index=["G1", "G2", "G3"]))
    path = _gene_set_session(tmp_path, ["G1"], data_mapping=files.mapping(expression=expression), var_count=var_count)
    with pytest.warns(UserWarning, match=re.escape(GENE_SETS_UNS_ONLY_NOTE)):
        report = read_session(path).apply(adata, dry_run=False)
    assert [row.status for row in report.plan.rows] == ["write"]
    assert adata.uns["cellucid"]["gene_sets"]["Markers"]["genes"] == ["G1"]
    assert list(adata.var.columns) == []


def test_a_gene_sets_genes_the_dataset_lacks_are_left_out_of_its_column(tmp_path: Path) -> None:
    adata = _small()
    path = _gene_set_session(tmp_path, ["G2", "NOPE", "X1", "X2", "X3"])
    message = 'The gene set "Markers" names 4 genes the dataset does not have (NOPE, X1, X2, …), left out of its var column.'
    with pytest.warns(UserWarning, match=re.escape(message)):
        read_session(path).apply(adata, dry_run=False)
    assert adata.var["cellucid_gene_set__Markers"].tolist() == [False, True]
    one = _small()
    with pytest.warns(UserWarning, match=re.escape('names 1 gene the dataset does not have (NOPE), left out')):
        read_session(_gene_set_session(tmp_path, ["NOPE"])).apply(one, dry_run=False)
    assert one.var["cellucid_gene_set__Markers"].tolist() == [False, False]


def test_a_gene_sets_column_reads_genes_as_the_data_mapping_names_them(tmp_path: Path) -> None:
    adata = _small()
    adata.var["symbol"] = ["CD3E", "LYZ"]
    path = _gene_set_session(tmp_path, ["LYZ"], data_mapping=files.mapping(geneNames="column:symbol"))
    read_session(path).apply(adata, dry_run=False)
    assert adata.var["cellucid_gene_set__Markers"].tolist() == [False, True]


def test_gene_set_columns_are_named_over_every_set_of_the_step() -> None:
    assert list(gene_set_columns(["T cells", "T-cells", "***", "T cells!"]).values()) == [
        "cellucid_gene_set__T_cells", "cellucid_gene_set__T_cells_2", "cellucid_gene_set__gene_set", "cellucid_gene_set__T_cells_3",
    ]


def test_a_gene_sets_var_column_reads_back_from_h5ad_and_zarr(golden: CellucidSession, tmp_path: Path) -> None:
    adata = golden_adata()
    golden.apply(adata, dry_run=False)
    expected = adata.var[T_CELL_MARKERS].tolist()
    adata.write_h5ad(tmp_path / "applied.h5ad")
    adata.write_zarr(tmp_path / "applied.zarr")
    for read in (ad.read_h5ad(tmp_path / "applied.h5ad"), ad.read_zarr(tmp_path / "applied.zarr")):
        assert read.var[T_CELL_MARKERS].dtype == bool and read.var[T_CELL_MARKERS].tolist() == expected
        # Read back, the same step is already there: nothing to write.
        assert not golden.apply(read, dry_run=False).applied


@pytest.mark.parametrize(
    "overlays",
    [
        pytest.param(files.overlays(user_fields=[files.user_category("user_cat_1", "cellucid_original_position", ["a"], [0] * 6)]),
                     id="user-field"),
        pytest.param(files.overlays(renamed_fields=[{"source": "obs", "key": "leiden", "name": "cellucid_original_position"}]),
                     id="renamed-onto"),
    ],
)
def test_an_apply_never_writes_or_renames_onto_the_priority_column(tmp_path: Path, overlays: dict[str, Any]) -> None:
    adata = _small()
    session = read_session(_write(tmp_path, {"set": {"field-overlays": overlays}}))
    with pytest.raises(ValueError, match=re.escape(ORIGINAL_POSITION_MESSAGE)):
        session.apply(adata)
    assert list(adata.obs.columns) == list(_small().obs.columns) and "cellucid" not in adata.uns


def test_an_apply_never_renames_or_relabels_the_priority_column(tmp_path: Path) -> None:
    # The column on an object not in a stored order (no priority record): it names no order, and still is never touched.
    for overlays in (
        files.overlays(renamed_fields=[{"source": "obs", "key": "cellucid_original_position", "name": "position"}]),
        files.overlays(renamed_categories=[{"source": "obs", "field": "cellucid_original_position", "category": "p0", "name": "first"}]),
    ):
        adata = _small()
        adata.obs["cellucid_original_position"] = pd.Categorical([f"p{i}" for i in range(adata.n_obs)])
        session = read_session(_write(tmp_path, {"set": {"field-overlays": overlays}}))
        with pytest.raises(ValueError, match=re.escape(ORIGINAL_POSITION_MESSAGE)), pytest.warns(UserWarning, match="names no order"):
            session.apply(adata)
        assert adata.obs["cellucid_original_position"].tolist() == [f"p{i}" for i in range(adata.n_obs)]


def test_presentation_domains_are_reported_not_applied(tmp_path: Path) -> None:
    camera = {"views": [{"zoom": 2.0, "ratio": 1.0}]}
    plan = read_session(_write(tmp_path, {"set": {"camera": camera}})).apply(_small())
    assert [(row.kind, row.domain, row.name, row.status, row.note) for row in plan.rows] == [
        ("presentation", "camera", "Camera", "not-applied", PRESENTATION_NOTE),
    ]
    report = plan.apply()
    assert not report.applied


def test_the_steps_expression_values_are_recorded_never_applied(tmp_path: Path) -> None:
    """The step's expression values go into its apply's record (what the person
    looked at the genes in); X keeps its values, and the domain is a presentation row."""
    gate = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [0, 1])))
    root = files.state(expression={"normalize": False, "log1p": False})
    shown = {"normalize": True, "log1p": True}
    path = _write(tmp_path, {"set": {"expression": shown}}, {"set": {"highlights": gate}}, root=root)
    adata = _small()
    stored = adata.X.copy()
    report = read_session(path).apply(adata, dry_run=False)
    assert report.applied
    assert [(row.kind, row.domain, row.name, row.status, row.note) for row in report.plan.rows if row.kind == "presentation"] == [
        ("presentation", "expression", "Expression values", "not-applied", PRESENTATION_NOTE),
    ]
    assert report.provenance is not None and report.provenance.expression == cellucid._history_model.ExpressionValues(normalize=True, log1p=True)
    assert [record["expression"] for record in history_applied(adata)] == [shown]
    assert np.array_equal(adata.X, stored)
    # Read back from a file, as a notebook reads it later.
    written = tmp_path / "applied.h5ad"
    adata.write_h5ad(written)
    assert cellucid.history_applied(ad.read_h5ad(written))[-1]["expression"] == shown
    # The step before (normalized and log1p asked, nothing written) records nothing.
    assert not read_session(path).apply(_small(), step="s1", dry_run=False).applied
    # A step whose expression values are the opening's records them too: values as stored.
    path = _write(tmp_path, {"set": {"highlights": gate}}, root=root)
    adata = _small()
    read_session(path).apply(adata, dry_run=False)
    assert history_applied(adata)[-1]["expression"] == {"normalize": False, "log1p": False}


def test_a_session_without_expression_values_records_none(tmp_path: Path) -> None:
    gate = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [0, 1])))
    adata = _small()
    report = read_session(_write(tmp_path, {"set": {"highlights": gate}})).apply(adata, dry_run=False)
    assert report.provenance is not None and report.provenance.expression is None
    assert history_applied(adata)[-1]["expression"] is None


@pytest.mark.parametrize(
    ("value", "message"),
    [
        ({"normalize": True}, r"expression value at step s1 has other fields than this version reads \(missing log1p\)"),
        ({"normalize": True, "log1p": 1}, r"expression value at step s1's log1p must be true or false"),
    ],
)
def test_expression_values_of_another_shape_are_refused(tmp_path: Path, value: dict[str, Any], message: str) -> None:
    gate = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [0, 1])))
    path = _write(tmp_path, {"set": {"expression": value, "highlights": gate}}, root=files.state(expression={"normalize": False, "log1p": False}))
    with pytest.raises(ValueError, match=message):
        read_session(path).apply(_small())


def test_a_steps_only_file_reads_and_is_not_applied(tmp_path: Path) -> None:
    path = _write(tmp_path, {"set": {"camera": {"views": [{"zoom": 2.0, "ratio": 1.0}]}}}, content="steps")
    session = read_session(path)
    assert (session.scope.content, session.scope.results) == ("steps", "recompute")
    assert session.tracker().steps["id"].tolist() == ["s0", "s1"]
    with pytest.raises(ValueError, match="shared as its steps only"):
        session.tracker().step("s1").state()
    with pytest.raises(ValueError, match="shared as its steps only, without the values they leave behind, so it cannot be applied"):
        session.apply(_small())


def test_a_state_only_file_reads_and_applies_as_any_file(tmp_path: Path) -> None:
    """Current state only: one tracker, the dataset as
    opened and one squashed step "Shared state" (no count, no author), with
    no note, turn, group or journal. It reads and applies as any file."""
    groups = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [0, 1])))
    shared = {"set": {"highlights": groups}, "label": "Shared state", "squashed": {"cause": "share", "count": None, "by": []}}
    session = read_session(_write(tmp_path, shared, content="state", times="exact"))
    assert (session.scope.content, session.scope.notes, session.scope.chats, session.scope.journal) == ("state", False, False, False)
    tracker = session.tracker()
    assert tracker.steps[["id", "kind", "type", "label"]].to_dict("records") == [
        {"id": "s0", "kind": "root", "type": None, "label": "Opened"},
        {"id": "s1", "kind": "squashed", "type": None, "label": "Shared state"},
    ]
    assert str(tracker.steps["v"].dtype) == "Int64" and tracker.steps["v"].isna().all(), "no action, no parameters version"
    assert (tracker.head_id, tracker.step("s1").squashed.cause, tracker.step("s1").squashed.by) == ("s1", "share", ())
    adata = _small()
    report = session.apply(adata, dry_run=False)
    assert report.applied and [(row.kind, row.name, row.cells) for row in report.plan.rows if row.kind == "highlight"] == [
        ("highlight", "cellucid_highlight__QC__Gate", 2),
    ]
    assert adata.obs["cellucid_highlight__QC__Gate"].tolist() == [True, True, False, False, False, False]


def test_a_state_file_holds_nothing_but_the_state(tmp_path: Path) -> None:
    groups = files.highlights(files.page("page_1", "QC", files.group("highlight_1", "Gate", [0, 1])))
    two_steps = [{"set": {"highlights": groups}}, {"set": {"camera": {"views": [{"zoom": 2.0, "ratio": 1.0}]}}}]
    with pytest.raises(ValueError, match=r"session\.scope of a state holds the root and one squashed step \(cause share\), shown"):
        read_session(_write(tmp_path, *two_steps, content="state"))
    # Privacy: a state shared without its history names nobody (no agent, model, notebook, person or turn).
    by_agent = {"set": {"highlights": groups}, "label": "Shared state", "agent": "Claude",
                "squashed": {"cause": "share", "count": None, "by": []}}
    with pytest.raises(ValueError, match=r"session\.scope of a state names nobody: every by is the person or the opening"):
        read_session(_write(tmp_path, by_agent, content="state"))


def test_the_options_are_exact(golden: CellucidSession) -> None:
    for option in ("dry_run", "overwrite", "keep_session"):
        with pytest.raises(TypeError, match=f"{option} must be exactly True or False"):
            apply_cellucid_session_to_anndata(golden_adata(), golden, **{option: 1})
    with pytest.raises(TypeError, match="step must be a string"):
        apply_cellucid_session_to_anndata(golden_adata(), golden, step=3)
    with pytest.raises(TypeError, match="session must be a CellucidSession"):
        apply_cellucid_session_to_anndata(golden_adata(), 3)


# =============================================================================
# What the share left out, and files that are damaged or newer
# =============================================================================


def _rewrite(path: Path, edit: Any, out: Path) -> Path:
    manifest, chunks = files.unframe(path.read_bytes())
    edit(manifest, chunks)
    out.write_bytes(files.frame(manifest, chunks))
    return out


def _chunk_index(manifest: dict[str, Any], chunk_id: str) -> int:
    return next(index for index, entry in enumerate(manifest["chunks"]) if entry["id"] == chunk_id)


def _edit_json_chunk(manifest: dict[str, Any], chunks: list[bytes], chunk_id: str, edit: Any) -> None:
    index = _chunk_index(manifest, chunk_id)
    entry = manifest["chunks"][index]
    payload = json.loads(files.payload(entry, chunks[index]))
    edit(payload)
    manifest["chunks"][index], chunks[index] = files.store(entry, json.dumps(payload, separators=(",", ":")).encode("utf-8"), entry["codec"])


def _flip_chunk_byte(manifest: dict[str, Any], chunks: list[bytes]) -> None:
    index = _chunk_index(manifest, "values/eager/0")
    body = bytearray(chunks[index])
    body[len(body) // 2] ^= 0xFF
    chunks[index] = bytes(body)


def _rename_page(manifest: dict[str, Any], chunks: list[bytes]) -> None:
    def edit(values: dict[str, Any]) -> None:
        for value in values.values():
            if isinstance(value, dict) and value.get("pages"):
                value["pages"][0]["name"] = "Page 9"

    _edit_json_chunk(manifest, chunks, "values/eager/0", edit)


def _extra_value(manifest: dict[str, Any], chunks: list[bytes]) -> None:
    value = {"unused": True}
    _edit_json_chunk(manifest, chunks, "values/lazy/0", lambda values: values.__setitem__(_hash_value(value), value))


def _drop_value(manifest: dict[str, Any], chunks: list[bytes]) -> None:
    _edit_json_chunk(manifest, chunks, "values/lazy/0", lambda values: values.pop(next(iter(values))))


def _extra_array(manifest: dict[str, Any], chunks: list[bytes]) -> None:
    array = files.typed("Uint8Array", [1, 2, 3])
    manifest["chunks"].insert(_chunk_index(manifest, "values/lazy/0"), {
        "id": f"array/{array.digest}", "kind": "array", "digest": array.digest, "encoding": "raw",
        "priority": "lazy", "codec": "gzip", "storedBytes": 0, "uncompressedBytes": 0,
    })
    chunks.insert(_chunk_index(manifest, f"array/{array.digest}"), files.stored(array.data.tobytes()))
    files.restore_sizes(manifest, chunks)


def _step_edit(edit: Any) -> Any:
    return lambda manifest, chunks: _edit_json_chunk(manifest, chunks, "steps/0/0", edit)


def _cut_gzip_stream(manifest: dict[str, Any], chunks: list[bytes]) -> None:
    index = _chunk_index(manifest, "values/eager/0")
    chunks[index] = chunks[index][: len(chunks[index]) // 2]
    manifest["chunks"][index]["storedBytes"] = len(chunks[index])


def _second_steps(steps: list[dict[str, Any]]) -> dict[str, Any]:
    return next(step for step in steps if step["id"] == "s1")


DAMAGES = [
    pytest.param(_cut_gzip_stream, 'the chunk "values/eager/0" is a gzip stream that ends early.', id="cut-gzip"),
    pytest.param(_flip_chunk_byte, 'the chunk "values/eager/0"', id="flipped-byte"),
    pytest.param(_rename_page, "its eager values holds a value that does not match its hash.", id="wrong-hash"),
    pytest.param(_extra_value, "it holds a value no step refers to.", id="unreferenced-value"),
    pytest.param(_drop_value, "a step refers to a value the file does not hold.", id="missing-value"),
    pytest.param(_extra_array, "it holds an array nothing refers to.", id="unreferenced-array"),
    pytest.param(lambda manifest, chunks: manifest["chunks"].append({**manifest["chunks"][-1], "kind": "picture", "hash": "a" * 64}),
                 "kind must be one of steps, parts, array, values, journal, chat.", id="picture-chunk"),
    pytest.param(_step_edit(lambda steps: _second_steps(steps).__setitem__("picture", None)), 'has an unknown field "picture"', id="picture-field"),
    pytest.param(_step_edit(lambda steps: _second_steps(steps).__setitem__("v", None)),
                 "records the version of its action's parameters (v)", id="no-version"),
    pytest.param(_step_edit(lambda steps: next(step for step in steps if step["id"] == "s3").__setitem__("group", None)),
                 'a group cuts the run of turn "py-1" between steps "s2" and "s3"; a group holds whole runs of a turn', id="split-turn"),
    pytest.param(_step_edit(lambda steps: next(step for step in steps if step["id"] == "s14").__setitem__("type", "color.by")),
                 "is squashed and records no action of its own", id="squashed-action"),
    pytest.param(lambda manifest, chunks: manifest["trackers"][1]["paused"].__setitem__("by", []),
                 "tracker.paused names who paused it (by)", id="paused-by-nobody"),
    pytest.param(lambda manifest, chunks: manifest["actions"].insert(0, manifest["actions"].pop(1)),
                 "session.actions must be sorted by type then version, each pair once", id="stamps-out-of-order"),
    pytest.param(lambda manifest, chunks: manifest["actions"][0].__setitem__("version", 2),
                 "session.actions must stamp exactly each action and version its steps use", id="stamp-unused"),
    pytest.param(lambda manifest, chunks: manifest["trackers"][1]["paused"].__setitem__("draft", {"camera": "a" * 64}),
                 "keeps a paused draft; a file squashes it into a step first", id="paused-draft"),
    pytest.param(lambda manifest, chunks: manifest["scope"].__setitem__("content", "state"),
                 "session.scope of a state holds one tracker of one branch", id="state-of-a-history"),
    pytest.param(lambda manifest, chunks: manifest["trackers"][0].__setitem__("headId", "s8"),
                 "head must be on its active branch's path.", id="head-off-branch"),
    pytest.param(lambda manifest, chunks: manifest["trackers"][0]["branches"][1].__setitem__("name", "MAIN"),
                 'branch "b0" and branch "b1" share the name "MAIN"', id="duplicate-branch-names"),
    pytest.param(lambda manifest, chunks: manifest["trackers"][0]["branches"][1].__setitem__("name", "B7"),
                 'A branch cannot be called "B7": that is the form of a branch id.', id="id-like-name"),
    pytest.param(_step_edit(lambda steps: _second_steps(steps)["children"].append("s5")), "lists a child that is not its own", id="tree"),
    pytest.param(_step_edit(lambda steps: _second_steps(steps)["touches"].append("camera")),
                 "must hold a value before and after for exactly the domains it touches", id="touches"),
    pytest.param(_step_edit(lambda steps: _second_steps(steps)["by"].__setitem__("origin", "robot")),
                 'by.origin must be one of "user", "agent", "host", "session", "replay", "startup", "outside"', id="origin"),
    pytest.param(_step_edit(lambda steps: _second_steps(steps).__setitem__("note", None)), 'has an unknown field "note"', id="old-step-key"),
    pytest.param(_step_edit(lambda steps: _second_steps(steps).__setitem__("branch", "b3")),
                 "steps must form one line starting right after its fork step", id="branch-membership"),
    pytest.param(lambda manifest, chunks: manifest.__setitem__("activeTracker", "t-0"), "its active tracker is not one it holds.", id="active"),
    pytest.param(lambda manifest, chunks: manifest.__setitem__("extra", 1), "must contain exactly", id="manifest-key"),
    pytest.param(lambda manifest, chunks: manifest["scope"].__setitem__("notes", False),
                 "holds notes or bookmarks its scope leaves out (notes).", id="notes-left-out"),
    pytest.param(lambda manifest, chunks: manifest["scope"].__setitem__("times", "dates"),
                 "holds a time of day at", id="dates-only"),
    pytest.param(lambda manifest, chunks: manifest["scope"].__setitem__("source", "name"),
                 "shared with the dataset's name only", id="source-name"),
    pytest.param(lambda manifest, chunks: manifest["scope"]["trackers"].reverse(),
                 "session.scope must name exactly the trackers it holds, in order", id="scope-trackers"),
    pytest.param(lambda manifest, chunks: manifest["chunks"][0].__setitem__("storedBytes", manifest["chunks"][0]["storedBytes"] + 1),
                 "is framed as", id="stored-bytes"),
    pytest.param(lambda manifest, chunks: manifest["chunks"].insert(0, manifest["chunks"].pop()),
                 'its chunk "steps/0/0" is out of order (eager before lazy; steps, parts, array, values, journal, chat).', id="order"),
    pytest.param(lambda manifest, chunks: manifest["requires"].__setitem__("obs", ["z", "a"]), "must be sorted and without duplicates.", id="requires"),
    pytest.param(lambda manifest, chunks: manifest.pop("chrome"), "missing: chrome", id="chrome-missing"),
    pytest.param(lambda manifest, chunks: manifest.__setitem__("chrome", []), "its page chrome must be exactly {pressed}.", id="chrome-not-record"),
    pytest.param(lambda manifest, chunks: manifest["chrome"].__setitem__("extra", 1), "its page chrome must be exactly {pressed}.", id="chrome-key"),
    pytest.param(lambda manifest, chunks: manifest["chrome"].__setitem__("pressed", []),
                 "its page chrome must map pressed-button groups to their pressed button.", id="chrome-pressed"),
    pytest.param(lambda manifest, chunks: manifest["chrome"]["pressed"].__setitem__("", "lasso"),
                 'its page chrome must name the pressed button of the group "" by its id.', id="chrome-empty-group"),
    pytest.param(lambda manifest, chunks: manifest["chrome"]["pressed"].__setitem__("highlight-mode", ""),
                 'its page chrome must name the pressed button of the group "highlight-mode" by its id.', id="chrome-empty-button"),
    pytest.param(lambda manifest, chunks: manifest["chrome"]["pressed"].__setitem__("highlight-mode", 3),
                 'its page chrome must name the pressed button of the group "highlight-mode" by its id.', id="chrome-not-text"),
]


@pytest.mark.core
@pytest.mark.parametrize(("damage", "message"), DAMAGES)
def test_a_damaged_file_is_refused_whole(tmp_path: Path, damage: Any, message: str) -> None:
    path = _rewrite(GOLDEN, damage, tmp_path / "damaged.cellucid-session")
    with pytest.raises(ValueError) as refused:
        read_session(path)
    assert str(refused.value).startswith("This session file is damaged: ")
    assert message in str(refused.value)


def test_a_file_from_a_newer_cellucid_is_newer_never_damaged(tmp_path: Path) -> None:
    def newer(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        manifest["app"]["version"] = "99.0.0"
        # A field this package's model does not know (a newer web app's).
        manifest["trackers"][0]["colour"] = "teal"

    path = _rewrite(GOLDEN, newer, tmp_path / "newer.cellucid-session")
    with pytest.raises(ValueError) as refused:
        read_session(path)
    assert str(refused.value) == newer_session_message("99.0.0", package_version())
    assert str(refused.value) == (
        f"This session was written by a newer Cellucid (99.0.0) than this package reads ({package_version()}); "
        "update the cellucid package to read it."
    )
    # A newer file this package reads is read.
    _rewrite(GOLDEN, lambda manifest, chunks: manifest["app"].__setitem__("version", "99.0.0"), tmp_path / "fine.cellucid-session")
    assert read_session(tmp_path / "fine.cellucid-session").app.version == "99.0.0"
    # The same unknown field from this version or an older one is damage.
    def older(manifest: dict[str, Any], chunks: list[bytes]) -> None:
        manifest["app"]["version"] = "0.9.0"
        manifest["trackers"][0]["colour"] = "teal"

    with pytest.raises(ValueError, match="^This session file is damaged: "):
        read_session(_rewrite(GOLDEN, older, tmp_path / "older.cellucid-session"))


def test_trailing_and_missing_bytes_are_refused(tmp_path: Path) -> None:
    data = GOLDEN.read_bytes()
    path = tmp_path / "framing.cellucid-session"
    path.write_bytes(data + b"\x00")
    with pytest.raises(ValueError, match="This session file is damaged: it holds bytes after its last chunk."):
        read_session(path)
    path.write_bytes(data[:-10])
    with pytest.raises(ValueError, match=r"This session file is damaged: the chunk .* only \d+ remain \(is the file truncated\?\)"):
        read_session(path)
    length = struct.unpack("<I", data[len(SESSION_SIGNATURE) : len(SESSION_SIGNATURE) + 4])[0]
    path.write_bytes(data[: len(SESSION_SIGNATURE) + 4 + length - 1])
    with pytest.raises(ValueError, match=f"its manifest length of {length} bytes is more than the file holds"):
        read_session(path)


#: The web app's own sentences (``bundle/format.js``), word for word.
WEB_EARLIER_FORMAT_3 = (
    "This session was saved by an earlier version of Cellucid (session format 3), which this version no "
    "longer reads: a session file now holds the history of every step (session format 4). It was not opened. "
    "Open it in the version of Cellucid that saved it, or re-create the session in this version."
)
WEB_NEWER_FORMAT_5 = (
    "This session was saved by a newer version of Cellucid (session format 5), which this version cannot "
    "read. Open it in that version of Cellucid."
)
WEB_NOT_A_SESSION = "This file is not a Cellucid session file of format 4 (it does not start with CELLUCID_SESSION/4)."


@pytest.mark.parametrize(
    ("prefix", "message"),
    [
        pytest.param(b"CELLUCID_SESSION/3\n", WEB_EARLIER_FORMAT_3, id="format-3"),
        pytest.param(b"CELLUCID_SESSION/5\n", WEB_NEWER_FORMAT_5, id="format-5"),
        pytest.param(b"CELLUCID_SESSION/0\n", WEB_NOT_A_SESSION, id="format-0"),
        pytest.param(b"PK\x03\x04 not a session at all", WEB_NOT_A_SESSION, id="other"),
    ],
)
def test_another_format_is_refused_with_the_web_apps_sentence(tmp_path: Path, prefix: bytes, message: str) -> None:
    path = tmp_path / "other.cellucid-session"
    path.write_bytes(prefix + b"\x10\x00\x00\x00{}" + b"\x00" * 64)
    with pytest.raises(ValueError) as refused:
        read_session(path)
    assert str(refused.value) == message
    with pytest.raises(ValueError) as refused:
        apply_cellucid_session_to_anndata(golden_adata(), path)
    assert str(refused.value) == message


def test_a_missing_file_is_not_found(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        read_session(tmp_path / "absent.cellucid-session")


def test_save_copies_the_file_whole(golden: CellucidSession, tmp_path: Path) -> None:
    out = golden.save(tmp_path / "copy.cellucid-session")
    assert out.read_bytes() == GOLDEN.read_bytes()
    assert golden.save(GOLDEN) == GOLDEN.resolve()


# =============================================================================
# Values hash exactly as the web app hashes them
# =============================================================================


@pytest.mark.parametrize(
    ("number", "text"),
    [
        (1e21, "1e+21"), (1e-7, "1e-7"), (123456789012345680000.0, "123456789012345680000"),
        (0.1 + 0.2, "0.30000000000000004"), (-1.5e-10, "-1.5e-10"), (5e-324, "5e-324"),
        (2**53 + 2, "9007199254740994"), (1 / 3, "0.3333333333333333"), (5.0, "5"), (100.0, "100"),
        (-0.001, "-0.001"), (1e-6, "0.000001"), (1.5e300, "1.5e+300"), (12345678901234567890, "12345678901234567000"),
    ],
)
def test_numbers_are_written_as_javascript_writes_them(number: float, text: str) -> None:
    assert _js_number(number) == text


def test_strings_and_keys_are_written_as_javascript_writes_them() -> None:
    assert _js_string(" \u007f\ud800\u0001\"\\") == '" \u007f\\ud800\\u0001\\"\\\\"'
    assert _js_stringify({"b": 1, "10": 2, "2": 3, "a": 4, "-1": 5, "01": 6}) == '{"2":3,"10":2,"b":1,"a":4,"-1":5,"01":6}'
    assert _hash_value({"\U0001f600": 1, "￿": 2}) == hashlib.sha256('{"\U0001f600":1,"￿":2}'.encode()).hexdigest()


def test_a_value_hashes_to_its_canonical_text() -> None:
    array = files.cells([1, 2])
    value = {"z": [1, 2.5, None, True, "t"], "a": array, "n": math.nan, "m": -0.0, "i": math.inf}
    text = (
        '{"a":{"$array":"' + array.digest + '"},"i":{"$number":"Infinity"},"m":{"$number":"-0"},'
        '"n":{"$number":"NaN"},"z":[1,2.5,null,true,"t"]}'
    )
    assert _hash_value(value) == hashlib.sha256(text.encode("utf-8")).hexdigest()
    assert _hash_value(copy.deepcopy({"b": 1, "a": 2})) == _hash_value({"a": 2, "b": 1})


# =============================================================================
# pull() of a session file
# =============================================================================


def test_pull_of_a_session_file_writes_the_state_it_shows(golden: CellucidSession) -> None:
    adata = golden_adata()
    summary = cellucid.pull(adata, GOLDEN)
    assert sorted(summary.columns) == sorted(["cellucid_highlight__highlight_1", "score", "cell_type (merged)"])
    # A pull writes results; renaming the object's own columns is apply's.
    assert adata.obs["cell_type"].cat.categories.tolist() == ["B cell", "T cell", "NK cell"]
    assert summary.gene_sets == {"T cell markers": ["CD3E", "NKG7"]}
    record = adata.uns["cellucid"]["pulled"]["cellucid_highlight__highlight_1"]
    assert (record["kind"], record["id"], record["pulled_from"]) == ("highlight", "highlight_1", GOLDEN.name)
    again = cellucid.pull(adata, GOLDEN)
    assert sorted(again.columns) == sorted(summary.columns)
    other = golden_adata()
    other.uns["cellucid"] = {"id": "another"}
    with pytest.raises(ValueError, match="saved on the dataset 'history_golden'"):
        cellucid.pull(other, GOLDEN)


# =============================================================================
# The format's closed sets are the web app's
# =============================================================================

def _web_source(relative: str) -> str:
    return web_file(f"assets/js/{relative}").read_text(encoding="utf-8")


def test_the_readers_constants_are_the_web_apps() -> None:
    """The chunk kinds and their keys are the session README's (the format-4
    contract), the typed arrays the value codec's; every record's keys are the
    history model's, which ``test_history_model.py`` compares with the web
    app's."""
    from cellucid import session as reader

    readme = _web_source("app/session/README.md")
    kinds = set(re.findall(r"^\| `[a-z]+/[^`]*` \| `([a-z]+)` \|", readme, re.M))
    assert kinds == set(reader._CHUNK_KEYS_BY_KIND)
    assert "`id`, `kind`, `priority`, `codec`, `storedBytes`,\n`uncompressedBytes`" in readme
    assert {"id", "kind", "priority", "codec", "storedBytes", "uncompressedBytes"} == reader._CHUNK_BASE_KEYS
    parts = _web_source("app/history/value-parts.js")
    types = re.search(r"ARRAY_TYPES = Object\.freeze\(Object\.fromEntries\(\[(.*?)\]\.map", parts, re.S)
    assert types is not None
    assert set(re.findall(r"\b([A-Z][A-Za-z0-9]*Array)\b", types.group(1))) == set(reader._ARRAY_TYPES)
    format_js = _web_source("app/session/bundle/format.js")
    assert f"export const SESSION_FORMAT_VERSION = {reader.SESSION_FORMAT_VERSION};" in format_js
