"""The shared data-mapping conformance fixture, run against the Python rules.

``tests/fixtures/data_mapping_conformance.json`` is kept byte-for-byte in
cellucid-python, cellucid-r and the web app, and each repository pins the same
SHA-256. The rules it pins -- making names unique, the best-guess ranking, the
embedding-column checks, canonical embedding ids, the default expression, the
vector-field key rule and the default vector bindings, the representations
enabled when a dataset opens, the "nothing to plot" text, and DataMapping
validation -- are the ones every reader and writer must apply identically, so a
change to any of them changes the fixture, the pin, and the other two
implementations together.

The Python package runs the sections whose rules it implements: making names
unique and the note that reports it, the embedding-column check, the
vector-field key rule (``prepare()``, the launch defaults and the transition
drift helper name and read vector fields with it), DataMapping validation, the
modality a multimodal dataset opens with, the check of a mapping's modality
against the dataset's, and the neighbour graph a dataset opens with
(``prepare()`` writes the one it chooses).
The other sections are rules of the web app.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from cellucid import data_mapping as dm

FIXTURE = Path(__file__).parent / "fixtures" / "data_mapping_conformance.json"

#: The same digest is pinned in cellucid-r (tests/testthat/test-data-mapping-
#: conformance.R) and in the web app. Update all three together.
FIXTURE_SHA256 = "6f7b7ff0e4392ea00c7e21c7a868ec68d3ec83a63558227f9f7d2c300a666c32"


def _fixture() -> dict[str, Any]:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def _cases(section: str) -> list[dict[str, Any]]:
    return _fixture()[section]


def test_the_fixture_is_the_pinned_shared_copy() -> None:
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256
    fixture = _fixture()
    assert fixture["version"] == 7
    assert list(fixture) == [
        "_about",
        "version",
        "makeUniqueNames",
        "uniqueNamesNote",
        "bestGuessTier",
        "embeddingSlots",
        "enabledEmbeddings",
        "embeddingColumns",
        "canonicalEmbeddingId",
        "defaultExpression",
        "vectorFieldKeys",
        "defaultVectorBindings",
        "nothingPlottable",
        "dataMapping",
        "defaultModality",
        "mappingModality",
        "defaultGraph",
    ]
    assert all(fixture[section] for section in list(fixture)[2:])


@pytest.mark.parametrize("case", _cases("makeUniqueNames"), ids=lambda case: case["name"])
def test_make_unique_names(case: dict[str, Any]) -> None:
    expected = case["expected"]
    if "error" in expected:
        with pytest.raises(dm.DataMappingError) as raised:
            dm.make_unique_names(case["names"], index=case["index"], join=case["join"])
        assert raised.value.code == expected["error"]
        return
    result = dm.make_unique_names(case["names"], index=case["index"], join=case["join"])
    assert result.names == expected["names"]
    assert [list(pair) for pair in result.renamed] == expected["renamed"]
    assert result.filled_from_index == expected["filledFromIndex"]
    assert len(set(result.names)) == len(result.names)


@pytest.mark.filterwarnings("ignore:Suffix used:UserWarning")
@pytest.mark.parametrize(
    "case",
    [case for case in _cases("makeUniqueNames") if "error" not in case["expected"]],
    ids=lambda case: case["name"],
)
def test_make_unique_names_is_anndata_make_index_unique(case: dict[str, Any]) -> None:
    """The rule is anndata's own, once empty names have taken their index name."""
    import anndata.utils as anndata_utils
    import pandas as pd
    filled = [
        index_name if name in ("", None) else name
        for name, index_name in zip(
            case["names"],
            case["index"] if case["index"] is not None else [None] * len(case["names"]),
            strict=True,
        )
    ]
    reference = anndata_utils.make_index_unique(pd.Index(filled, dtype=object), join=case["join"])
    assert list(reference) == case["expected"]["names"]


@pytest.mark.parametrize("case", _cases("uniqueNamesNote"))
def test_unique_names_note(case: dict[str, Any]) -> None:
    assert (
        dm.unique_names_note(
            filled_from_index=case["filledFromIndex"],
            renamed=case["renamed"],
            source=case["source"],
        )
        == case["expected"]
    )


@pytest.mark.parametrize("case", _cases("embeddingColumns"))
def test_embedding_columns(case: dict[str, Any]) -> None:
    expected = case["expected"]
    if "error" in expected:
        with pytest.raises(dm.DataMappingError) as raised:
            dm.validate_embedding_columns(
                case["columns"], dimension=case["dimension"], width=case["width"]
            )
        assert raised.value.code == expected["error"]
        return
    assert (
        dm.validate_embedding_columns(
            case["columns"], dimension=case["dimension"], width=case["width"]
        )
        == expected["columns"]
    )


@pytest.mark.parametrize("case", _cases("vectorFieldKeys"), ids=lambda case: case["name"])
def test_vector_field_keys(case: dict[str, Any]) -> None:
    found = dm.vector_field_keys(case["keys"])
    assert [
        {"key": key, "field": parsed.field, "basis": parsed.basis, "embeddingKey": parsed.embedding_key}
        for key, parsed in found.items()
    ] == case["expected"]
    keys = set(case["keys"])
    for key in case["keys"]:
        parsed = dm.parse_vector_field_key(key, keys)
        if key in found:
            assert parsed == found[key]
            assert parsed is not None
            assert f"{parsed.field}_{parsed.basis}" == key
            # The embedding is in the field's own namespace (a modality's ``rna:``).
            assert parsed.embedding_key == dm.obsm_namespace(parsed.field)[0] + "X_" + parsed.basis
        else:
            assert parsed is None


def test_vector_field_key_rule_splits_at_any_character() -> None:
    keys = {"X_umap", "X_\U0001F600"}
    assert dm.parse_vector_field_key("a\r_umap", keys) == dm.VectorFieldKey("a\r", "umap")
    assert dm.parse_vector_field_key("a\n b_umap", keys) == dm.VectorFieldKey("a\n b", "umap")
    assert dm.parse_vector_field_key("\U0001F600_\U0001F600", keys) == dm.VectorFieldKey(
        "\U0001F600", "\U0001F600"
    )
    assert dm.parse_vector_field_key("velocity_umap\n", keys) is None
    assert dm.parse_vector_field_key("X_umap", keys) is None
    assert dm.parse_vector_field_key(3, keys) is None


@pytest.mark.parametrize("case", _cases("dataMapping"), ids=lambda case: case["name"])
def test_data_mapping_validation(case: dict[str, Any]) -> None:
    expected = case["expected"]
    if "error" in expected:
        with pytest.raises(dm.DataMappingError) as raised:
            dm.validate_data_mapping(case["value"])
        assert raised.value.code == expected["error"]
        return
    canonical = dm.validate_data_mapping(case["value"])
    assert canonical == expected
    # The canonical form is a fixed point and is serialized key-for-key.
    assert dm.validate_data_mapping(canonical) == canonical
    assert list(canonical) == list(expected)


@pytest.mark.parametrize("case", _cases("defaultGraph"), ids=lambda case: case["name"])
def test_default_graph(case: dict[str, Any]) -> None:
    assert dm.default_graph(case["keys"], modality=case["modality"]) == case["expected"]


@pytest.mark.parametrize("case", _cases("defaultModality"))
def test_default_modality(case: dict[str, Any]) -> None:
    assert dm.default_modality(case["names"]) == case["expected"]


@pytest.mark.parametrize("case", _cases("mappingModality"), ids=lambda case: case["name"])
def test_mapping_modality(case: dict[str, Any]) -> None:
    if case["expected"] is None:
        dm.require_mapping_modality(case["mapping"], case["modalities"])
        return
    with pytest.raises(dm.DataMappingError) as raised:
        dm.require_mapping_modality(case["mapping"], case["modalities"])
    assert raised.value.code == case["expected"]["error"]


def test_the_embeddings_shown_at_open_follow_the_web_rule_from_metadata() -> None:
    """``enabled_embedding_keys`` is the web app's ``resolveEmbeddingSlots`` and
    ``defaultEnabledEmbeddings`` without its load checks."""

    def entries(**widths: int | None) -> list[dm.ObsmEntry]:
        return [
            dm.ObsmEntry(key=key, width=width, usable=width is not None)
            for key, width in widths.items()
        ]

    # The default rule fills 2D; narrow arrays are shown beside it, a wide one
    # and vector keys are not.
    assert dm.enabled_embedding_keys(
        entries(X_umap=2, X_tsne=3, X_pca=50, velocity_umap=2),
        named_slots={},
        latent_key=None,
    ) == {"X_umap", "X_tsne"}
    # The latent space is never shown beside.
    assert dm.enabled_embedding_keys(
        entries(X_umap=2, X_tsne=2), named_slots={}, latent_key="X_tsne"
    ) == {"X_umap"}
    # Nothing narrow: the best guess shows the first wide array, by tier.
    assert dm.enabled_embedding_keys(
        entries(X_pca=50, X_draw_graph_fa=8), named_slots={}, latent_key=None
    ) == {"X_draw_graph_fa"}
    # A named slot is shown whatever its width.
    assert dm.enabled_embedding_keys(
        entries(X_umap=2, X_pca=50), named_slots={"2d": "X_pca"}, latent_key=None
    ) == {"X_umap", "X_pca"}
    assert dm.best_guess_order(["X_pca", "X_spatial", "X_other", "X_tsne"]) == [
        "X_tsne",
        "X_spatial",
        "X_other",
        "X_pca",
    ]
