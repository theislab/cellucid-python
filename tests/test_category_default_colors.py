"""Cellucid's own category colors, as the web app gives them.

A user field a session applies without colors of its own is written with the
colors the viewer showed it in: Cellucid's category colors by position
(``data/palettes.js`` ``getCategoryColor``). The shared fixture, written by
the web app's ``scripts/generate-category-default-colors.mjs``, pins the
first 1,024 byte for byte in cellucid, cellucid-python and cellucid-r.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from cellucid._category_colors import default_category_color

FIXTURE = Path(__file__).parent / "fixtures" / "category_default_colors.json"
FIXTURE_SHA256 = "66e652c560afa8aa5bbe149ed87379fd9b937e5a3ca0585d669fdc09092576e0"


def test_the_colors_are_the_web_apps() -> None:
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256
    expected = json.loads(FIXTURE.read_text(encoding="utf-8"))["colors"]
    assert [default_category_color(index) for index in range(len(expected))] == expected


@pytest.mark.parametrize("index", [-1, 1.0, True, "0"])
def test_a_place_is_a_non_negative_whole_number(index: object) -> None:
    with pytest.raises(ValueError, match="non-negative whole number"):
        default_category_color(index)  # type: ignore[arg-type]
