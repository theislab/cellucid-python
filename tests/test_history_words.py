"""A history in words: the methods paragraph and the pseudocode, by the web app's rules.

The golden history (``tests/fixtures/history_golden``) pins the text of three
branches of its annotation tracker byte for byte (``test_sessions.py``). These
tests cover the wording rules the golden branches do not reach, on the golden
tracker's own steps (the Scratch branch up to its paused period, s14) changed
with :func:`dataclasses.replace`: a squashed step of each cause (recording
paused, steps merged, the history cleared, a state shared) listed in the
methods and kept as a pseudocode comment, the sentence that says a squashed
step holds only its result, everyone who acted while recording was paused,
and the "Before"/"After" headings around a group. The same strings are
asserted by R's ``test-history-words.R``.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from cellucid import _history_model as model
from cellucid.history import methods_paragraph, pseudocode_text, squash_phrase
from cellucid.session import read_session

GOLDEN = Path(__file__).parent / "fixtures" / "history_golden" / "history_golden.cellucid-session"
DATASET = {"name": "History golden", "cells": 600, "genes": 4}
OPENING = (
    "Single-cell data were explored interactively in Cellucid. The dataset “History golden” (600 cells and 4 genes) was"
    " opened in full."
)
COMPLETE = "The complete history, with every parameter, can be reopened from the accompanying Cellucid session file."


@pytest.fixture(scope="module")
def tracker():  # type: ignore[no-untyped-def]
    return read_session(GOLDEN).tracker("Annotation")


def test_squash_phrases() -> None:
    def phrase(cause: str, count: int | None) -> str:
        return squash_phrase(model.Squash(cause=cause, count=count, by=(), carried=model.Carried(notes=(), bookmarks=())))

    assert phrase("paused", None) == "Changes made with recording paused"
    assert phrase("steps", 3) == "3 steps merged into one"
    assert phrase("clear", 9) == "A fresh start from the state reached before the history was cleared"
    assert phrase("share", None) == "A shared state, received without its history"


def test_a_paused_period_that_changed_the_data_is_listed_and_said_plainly(tracker) -> None:  # type: ignore[no-untyped-def]
    scratch = list(tracker.path_to("s14"))
    # Recording was paused while the filters and the colors changed (the person and a notebook acted).
    scratch[-1] = dataclasses.replace(scratch[-1], touches=("filters", "colors"))
    text = methods_paragraph(scratch, dataset=DATASET, version=None, notes=tracker.notes[:0], groups=tracker.groups)
    assert text == "\n\n".join([
        f"{OPENING} The 6 steps applied are listed below, in order."
        " For step 6, only the state it reached was kept, not the operations that led to it."
        " 1 further step only changed the camera, the panel layout, colors or display settings and is not listed."
        f" {COMPLETE}",
        "Before QC\nThe following 1 step was applied: (1) renamed NK cell to NK in cell_type.",
        'QC\nThe following 2 steps were applied, in this order: (2) highlighted 40 cells as "Doublets" on page "QC";'
        " (3) added the continuous column score. Steps 2–3 were taken by a Python notebook through the Cellucid viewer API.",
        "After QC\nThe following 3 steps were applied, in this order: (4) ran Marker genes; (5) kept the gene set T cell"
        " markers (2 genes); (6) changes made with recording paused (filters, colors). Steps 4–5 were taken by an AI"
        " assistant (Claude, claude-sonnet-4-6) at the user's request, and step 6 by a Python notebook through the"
        " Cellucid viewer API.",
    ])


def test_each_squash_cause_in_the_methods_and_the_pseudocode(tracker) -> None:  # type: ignore[no-untyped-def]
    scratch = list(tracker.path_to("s14"))
    root, paused = scratch[0], scratch[-1]
    person = paused.by
    cases = {
        "steps": (3, "3 steps merged into one"),
        "clear": (5, "a fresh start from the state reached before the history was cleared"),
        "share": (None, "a shared state, received without its history"),
    }
    for cause, (count, said) in cases.items():
        step = dataclasses.replace(
            paused, parent=root.id, group=None, bookmark=None, touches=("field-overlays", "highlights"),
            squashed=model.Squash(
                cause=cause, count=count, by=() if cause == "share" else (person,), carried=model.Carried(notes=(), bookmarks=()),
            ),
        )
        assert methods_paragraph([root, step], dataset=DATASET, version=None) == (
            f"{OPENING} The following 1 step was applied: (1) {said} (fields and labels, highlights)."
            f" For step 1, only the state it reached was kept, not the operations that led to it. {COMPLETE}"
        )
        assert pseudocode_text([root, step], dataset_name=None, version=None, tracker="T", branch="B") == (
            '# Cellucid\n# Tracker "T", branch "B": 1 step\n'
            f"# {said[0].upper()}{said[1:]}: Fields and labels, Highlights (values only)\n"
        )


def test_a_squashed_step_that_only_presented_is_counted(tracker) -> None:  # type: ignore[no-untyped-def]
    path = list(tracker.path_to("s14"))
    text = methods_paragraph(path, dataset=DATASET, version=None, groups=tracker.groups)
    assert "2 further steps only changed the camera, the panel layout, colors or display settings and are not listed." in text
    assert "recording paused" not in text
    assert pseudocode_text(path, dataset_name=None, version=None, tracker="T", branch="B", groups=tracker.groups).endswith(
        "# Changes made with recording paused: Colors, Settings (values only)\n"
    )


@pytest.mark.parametrize(
    ("open_cells", "share"),
    [(2_000, "10%"), (1_100, "5.5%"), (19_999, "99%"), (3, "0.015%"), (19_900, "99%")],
)
def test_a_share_of_the_cells_is_written_as_the_web_app_writes_it(tracker, open_cells: int, share: str) -> None:  # type: ignore[no-untyped-def]
    # formatShare (cellucid/assets/js/utils/format-units.js): whole from 10%, never 100% short of
    # every cell; two significant figures below. Before round 9 this said "6%" for 5.5% and the
    # web app "5.5%" (R9U-3), and below 1% the web app's own two copies differed.
    dataset = {**DATASET, "cells": 20_000, "open_cells": {"open": open_cells, "of": 20_000, "method": "uniform", "seed": 0}}
    text = methods_paragraph(list(tracker.path_to("s14"))[:1], dataset=dataset, version=None)
    assert f"was opened on {open_cells:,} of its 20,000 cells ({share})" in text
