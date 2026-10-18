"""The session fingerprint requires a cell-order record naming its embedding.

The viewer records which cell ordering a session was saved against, so a dataset
republished with the same name and counts but a different row order cannot
silently re-point every saved selection at different cells. The digest is over
the coordinates of one embedding, so the record names that embedding: an
exported embedding id with ``columns: null``, or an AnnData ``obsm`` key with
the 0-based columns read. The record and the data mapping are both required;
the viewer refuses a session without them and so does this reader. The digest
itself cannot be re-derived here, which is the viewer's job; its shape is.
"""

from __future__ import annotations

import copy
from typing import Any

import pytest

from cellucid.session_bundle import _validate_fingerprint

_EXPORT_MAPPING: dict[str, Any] = {
    "kind": "export",
    "embeddings": {
        "1d": None,
        "2d": {"key": "X_tsne", "columns": None},
        "3d": {"key": "X_umap_3d", "columns": None},
    },
    "geneNames": "index",
    "vectorFields": [],
}

_ANNDATA_MAPPING: dict[str, Any] = {
    "kind": "anndata",
    "embeddings": {"1d": None, "2d": {"key": "X_pca", "columns": [3, 1]}, "3d": None},
    "knn": None,
    "geneNames": "column:gene_symbols",
    "expression": "X",
    "latent": None,
    "vectorFields": [],
}


def _fingerprint(*, mapping: dict[str, Any] | None = None, **overrides: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "sourceType": "local-demo",
        "datasetId": "pancreas",
        "cellCount": 3696,
        "varCount": 3753,
        "cellOrder": {
            "embedding": {"key": "X_tsne", "columns": None},
            "dimension": 2,
            "digest": "0123456789abcdef",
        },
        "dataMapping": copy.deepcopy(_EXPORT_MAPPING if mapping is None else mapping),
    }
    base.update(overrides)
    return base


def test_a_complete_export_fingerprint_is_accepted() -> None:
    validated = _validate_fingerprint(_fingerprint())
    assert validated["cellOrder"]["embedding"] == {"key": "X_tsne", "columns": None}
    assert validated["dataMapping"] == _EXPORT_MAPPING


def test_a_complete_anndata_fingerprint_is_accepted() -> None:
    fingerprint = _fingerprint(
        mapping=_ANNDATA_MAPPING,
        cellOrder={
            "embedding": {"key": "X_pca", "columns": [3, 1]},
            "dimension": 2,
            "digest": "fedcba9876543210",
        },
    )
    assert _validate_fingerprint(fingerprint)["dataMapping"] == _ANNDATA_MAPPING


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_every_embedding_dimension_is_accepted(dimension: int) -> None:
    _validate_fingerprint(
        _fingerprint(
            cellOrder={
                "embedding": {"key": "X_tsne", "columns": None},
                "dimension": dimension,
                "digest": "a" * 16,
            }
        )
    )


def test_the_cell_order_embedding_need_not_be_a_slot_default() -> None:
    # A session saved while a view showed another layout digests that layout.
    _validate_fingerprint(
        _fingerprint(
            cellOrder={
                "embedding": {"key": "X_umap", "columns": None},
                "dimension": 2,
                "digest": "a" * 16,
            }
        )
    )


_EMBEDDING = {"key": "X_tsne", "columns": None}


@pytest.mark.parametrize(
    ("cell_order", "message"),
    [
        ("not-an-object", "must be an object"),
        ({"embedding": _EMBEDDING, "dimension": 2}, "missing digest"),
        ({"embedding": _EMBEDDING, "digest": "a" * 16}, "missing dimension"),
        ({"dimension": 2, "digest": "a" * 16}, "missing embedding"),
        (
            {"embedding": _EMBEDDING, "dimension": 2, "digest": "a" * 16, "extra": 1},
            "unknown extra",
        ),
        ({"embedding": _EMBEDDING, "dimension": 0, "digest": "a" * 16}, "must be exactly 1, 2, or 3"),
        ({"embedding": _EMBEDDING, "dimension": 4, "digest": "a" * 16}, "must be exactly 1, 2, or 3"),
        (
            {"embedding": _EMBEDDING, "dimension": True, "digest": "a" * 16},
            "must be exactly 1, 2, or 3",
        ),
        ({"embedding": _EMBEDDING, "dimension": 2, "digest": "A" * 16}, "16 lowercase hexadecimal"),
        ({"embedding": _EMBEDDING, "dimension": 2, "digest": "a" * 15}, "16 lowercase hexadecimal"),
        ({"embedding": _EMBEDDING, "dimension": 2, "digest": "a" * 17}, "16 lowercase hexadecimal"),
        ({"embedding": _EMBEDDING, "dimension": 2, "digest": "g" * 16}, "16 lowercase hexadecimal"),
        ({"embedding": _EMBEDDING, "dimension": 2, "digest": 1}, "16 lowercase hexadecimal"),
        (
            {"embedding": "X_tsne", "dimension": 2, "digest": "a" * 16},
            "cellOrder.embedding must be an object",
        ),
        (
            {"embedding": {"key": "X_tsne"}, "dimension": 2, "digest": "a" * 16},
            "cellOrder.embedding fields: missing columns",
        ),
        (
            {"embedding": {"key": "", "columns": None}, "dimension": 2, "digest": "a" * 16},
            "cellOrder.embedding.key must be a non-empty string",
        ),
        (
            {"embedding": {"key": "X_tsne", "columns": [0, 1]}, "dimension": 2, "digest": "a" * 16},
            "cellOrder.embedding.columns must be null",
        ),
    ],
)
def test_a_malformed_export_cell_order_is_rejected(cell_order: Any, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _validate_fingerprint(_fingerprint(cellOrder=cell_order))


@pytest.mark.parametrize(
    ("embedding", "dimension", "message"),
    [
        ({"key": "X_pca", "columns": None}, 2, "cellOrder.embedding.columns must be a list"),
        ({"key": "X_pca", "columns": [0, 1, 2]}, 2, "must list exactly 2 columns"),
        ({"key": "X_pca", "columns": [1, 1]}, 2, "each axis needs a different column"),
        ({"key": "X_pca", "columns": [0, -1]}, 2, "non-negative integers"),
        ({"key": "X_pca", "columns": [0, True]}, 2, "non-negative integers"),
    ],
)
def test_an_anndata_cell_order_lists_one_column_per_axis(
    embedding: dict[str, Any], dimension: int, message: str
) -> None:
    fingerprint = _fingerprint(
        mapping=_ANNDATA_MAPPING,
        cellOrder={"embedding": embedding, "dimension": dimension, "digest": "a" * 16},
    )
    with pytest.raises(ValueError, match=message):
        _validate_fingerprint(fingerprint)


@pytest.mark.parametrize("missing", ["cellOrder", "dataMapping", "varCount"])
def test_every_fingerprint_field_is_required(missing: str) -> None:
    incomplete = _fingerprint()
    del incomplete[missing]
    with pytest.raises(ValueError, match=f"missing {missing}"):
        _validate_fingerprint(incomplete)


def test_an_unknown_fingerprint_field_is_still_rejected() -> None:
    with pytest.raises(ValueError, match="unknown somethingElse"):
        _validate_fingerprint(_fingerprint(somethingElse=1))


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"kind": "h5ad"}, "dataMapping.kind must be"),
        ({"knn": None}, "dataMapping fields: unknown knn"),
        ({"geneNames": "symbols"}, "dataMapping.geneNames must be"),
        (
            {"embeddings": {"1d": None, "2d": None, "3d": None}},
            "dataMapping.embeddings must fill at least one slot",
        ),
        (
            {"vectorFields": [{"field": "velocity", "embedding": {"key": "X_tsne", "columns": None}}] * 2},
            "dataMapping.vectorFields repeats a binding",
        ),
    ],
)
def test_a_malformed_data_mapping_is_rejected(change: dict[str, Any], message: str) -> None:
    mapping = copy.deepcopy(_EXPORT_MAPPING)
    mapping.update(change)
    with pytest.raises(ValueError, match=message):
        _validate_fingerprint(_fingerprint(mapping=mapping))


def test_the_returned_fingerprint_is_a_copy() -> None:
    fingerprint = _fingerprint()
    validated = _validate_fingerprint(fingerprint)
    validated["dataMapping"]["embeddings"]["2d"]["key"] = "changed"
    assert fingerprint["dataMapping"]["embeddings"]["2d"]["key"] == "X_tsne"
