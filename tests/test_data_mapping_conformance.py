"""The shared data-mapping conformance fixture, run against the Python rules.

``tests/fixtures/data_mapping_conformance.json`` is kept byte-for-byte in
cellucid-python, cellucid-r and the web app, and each repository pins the same
SHA-256. The rules it pins -- making names unique, the best-guess ranking, the
embedding-column checks, canonical embedding ids, the default expression and
vector bindings, the "nothing to plot" text, and DataMapping validation -- are
the ones every reader and writer must apply identically, so a change to any of
them changes the fixture, the pin, and the other two implementations together.
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
FIXTURE_SHA256 = "da5f48332d4370bc74ad162e129e854398b1392ad76f7bdf68b389ad0bd11b65"


def _fixture() -> dict[str, Any]:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def _cases(section: str) -> list[dict[str, Any]]:
    return _fixture()[section]


def test_the_fixture_is_the_pinned_shared_copy() -> None:
    assert hashlib.sha256(FIXTURE.read_bytes()).hexdigest() == FIXTURE_SHA256
    fixture = _fixture()
    assert fixture["version"] == 1
    assert list(fixture) == [
        "_about",
        "version",
        "makeUniqueNames",
        "uniqueNamesNote",
        "bestGuessTier",
        "embeddingSlots",
        "embeddingColumns",
        "canonicalEmbeddingId",
        "defaultExpression",
        "defaultVectorBindings",
        "nothingPlottable",
        "dataMapping",
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
    pd = pytest.importorskip("pandas")
    anndata_utils = pytest.importorskip("anndata.utils")
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


@pytest.mark.parametrize("case", _cases("bestGuessTier"), ids=lambda case: case["key"])
def test_best_guess_tier(case: dict[str, Any]) -> None:
    assert dm.best_guess_tier(case["key"]) == case["tier"]


@pytest.mark.parametrize("case", _cases("embeddingSlots"), ids=lambda case: case["name"])
def test_embedding_slots(case: dict[str, Any]) -> None:
    entries = [
        dm.ObsmEntry(key=entry["key"], width=entry["width"], usable=entry["usable"])
        for entry in case["obsm"]
    ]
    failures = case["loadFailures"]
    checked: list[tuple[str, list[int]]] = []

    def load_check(key: str, columns: list[int]) -> str | None:
        checked.append((key, list(columns)))
        return failures.get(key)

    guess = dm.resolve_embedding_slots(
        entries,
        latent_key=case["latentKey"],
        load_check=load_check,
    )
    assert guess.slots == case["expected"]["slots"]
    assert guess.guessed == case["expected"]["guessed"]
    assert [list(item) for item in guess.rejected] == case["expected"]["rejected"]
    # A candidate is load-checked at most once, so the guess is bounded by the
    # candidate count.
    assert len({key for key, _ in checked}) == len(checked)


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


@pytest.mark.parametrize("case", _cases("canonicalEmbeddingId"))
def test_canonical_embedding_id(case: dict[str, Any]) -> None:
    assert (
        dm.canonical_embedding_id(case["key"], case["columns"], width=case["width"])
        == case["expected"]
    )


@pytest.mark.parametrize("case", _cases("defaultExpression"))
def test_default_expression(case: dict[str, Any]) -> None:
    value, guessed = dm.default_expression(
        has_x=case["hasX"],
        has_raw=case["hasRaw"],
        layers=case["layers"],
    )
    assert {"expression": value, "guessed": guessed} == case["expected"]


@pytest.mark.parametrize("case", _cases("defaultVectorBindings"))
def test_default_vector_bindings(case: dict[str, Any]) -> None:
    assert dm.default_vector_bindings(case["keys"], case["slots"]) == case["expected"]


@pytest.mark.parametrize("case", _cases("nothingPlottable"))
def test_nothing_plottable_message(case: dict[str, Any]) -> None:
    descriptions = [dm.describe_obsm_entry(**entry) for entry in case["obsm"]]
    assert dm.nothing_plottable_message(n_obs=case["nObs"], descriptions=descriptions) == case[
        "expected"
    ]


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


def test_the_best_guess_skips_vector_keys_and_the_latent_space_by_role() -> None:
    entries = [
        dm.ObsmEntry(key="velocity_umap_2d", width=2, usable=True),
        dm.ObsmEntry(key="X_scvi", width=2, usable=True),
    ]
    guess = dm.resolve_embedding_slots(entries, latent_key="X_scvi")
    assert guess.slots == {"1d": None, "2d": None, "3d": None}
    described = [
        {"key": "velocity_umap_2d", "storage": "dense", "shape": [5, 2], "usable": True, "reason": None},
        {"key": "X_scvi", "storage": "dense", "shape": [5, 2], "usable": True, "reason": None},
    ]
    defects = [
        dm.obsm_entry_defect(item, latent_key="X_scvi", rejected={}) for item in described
    ]
    assert defects == ["a vector field", "the latent space"]
