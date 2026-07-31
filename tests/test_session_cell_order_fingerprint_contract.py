"""The session reader accepts a fingerprint with or without ``cellOrder``.

The viewer records which cell ordering a session was saved against, so a dataset
republished with the same name and counts but a different row order cannot
silently re-point every saved selection at different cells. This reader accepts
both shapes: it cannot verify the digest -- which is taken over exported
coordinates -- so enforcing that identity is the viewer's job, and accepting
both is what lets the viewer add the field without stranding existing files.
"""

from __future__ import annotations

from typing import Any

import pytest

from cellucid.session_bundle import _validate_fingerprint


def _fingerprint(**overrides: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "sourceType": "local-demo",
        "datasetId": "pancreas",
        "cellCount": 3696,
        "varCount": 3753,
    }
    base.update(overrides)
    return base


def test_a_fingerprint_without_cell_order_is_accepted() -> None:
    _validate_fingerprint(_fingerprint())


def test_a_fingerprint_with_cell_order_is_accepted() -> None:
    _validate_fingerprint(_fingerprint(cellOrder={"dimension": 2, "digest": "0123456789abcdef"}))


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_every_embedding_dimension_is_accepted(dimension: int) -> None:
    _validate_fingerprint(_fingerprint(cellOrder={"dimension": dimension, "digest": "a" * 16}))


@pytest.mark.parametrize(
    ("cell_order", "message"),
    [
        ("not-an-object", "must be an object"),
        ({"dimension": 2}, "missing digest"),
        ({"digest": "a" * 16}, "missing dimension"),
        ({"dimension": 2, "digest": "a" * 16, "extra": 1}, "unknown extra"),
        ({"dimension": 0, "digest": "a" * 16}, "must be exactly 1, 2, or 3"),
        ({"dimension": 4, "digest": "a" * 16}, "must be exactly 1, 2, or 3"),
        ({"dimension": True, "digest": "a" * 16}, "must be exactly 1, 2, or 3"),
        ({"dimension": 2, "digest": "A" * 16}, "16 lowercase hexadecimal"),
        ({"dimension": 2, "digest": "a" * 15}, "16 lowercase hexadecimal"),
        ({"dimension": 2, "digest": "a" * 17}, "16 lowercase hexadecimal"),
        ({"dimension": 2, "digest": "g" * 16}, "16 lowercase hexadecimal"),
        ({"dimension": 2, "digest": 1}, "16 lowercase hexadecimal"),
    ],
)
def test_a_malformed_cell_order_is_rejected(cell_order: Any, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _validate_fingerprint(_fingerprint(cellOrder=cell_order))


def test_an_unknown_fingerprint_field_is_still_rejected() -> None:
    with pytest.raises(ValueError, match="unknown somethingElse"):
        _validate_fingerprint(_fingerprint(somethingElse=1))


def test_a_missing_required_field_is_still_rejected() -> None:
    incomplete = _fingerprint()
    del incomplete["varCount"]
    with pytest.raises(ValueError, match="missing varCount"):
        _validate_fingerprint(incomplete)
