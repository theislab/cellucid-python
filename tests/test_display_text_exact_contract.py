"""Every user-supplied string the viewer prints must read as the value it stores.

A dataset name, a description, a citation and a source line are drawn as
stored in the viewer. A character with no glyph inside them changes the stored
value without changing the picture, so ``prepare()`` rejects such text instead
of removing it. White space at either end of these descriptive texts is the
one exception: it is removed, as the web app reads it; a URL keeps the
strict rule.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cellucid import prepare

_EMBEDDING = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32)


def _prepare_kwargs():
    return {
        "embeddings": {"X_umap": _EMBEDDING.copy()},
        "dataset_id": "display-text",
        "dataset_name": "Display text",
    }


_PROSE = ["dataset_name", "dataset_description", "citation_text", "source_name"]
_URLS = ["citation_url", "source_url"]
#: Refused in every displayed text: a character with no glyph inside it.
_INSIDE = [
    ("hid\u200bden", "contains U+200B ZERO WIDTH SPACE"),
    ("cont\x01rol", "contains U+0001 (control character)"),
    ("Li\ud800ver", "contains U+D800 (unpaired surrogate)"),
]
#: Refused at the ends of a URL; removed from the ends of descriptive text.
_EDGES = [
    (" padded ", "starts with U+0020 SPACE"),
    ("padded\u00a0", "ends with U+00A0 NO-BREAK SPACE"),
    ("\ufeffLiver", "contains U+FEFF ZERO WIDTH NO-BREAK SPACE"),
]


@pytest.mark.parametrize(
    ("argument", "value", "expected"),
    [(argument, value, expected) for argument in [*_PROSE, *_URLS] for value, expected in _INSIDE]
    + [(argument, value, expected) for argument in _URLS for value, expected in _EDGES],
)
def test_prepare_rejects_undisplayable_identity_text_before_output_mutation(
    tmp_path,
    argument,
    value,
    expected,
):
    out_path = tmp_path / "must-not-exist.h5ad"
    kwargs = _prepare_kwargs()
    if argument == "citation_url":
        kwargs["citation_text"] = "Exact citation"
    if argument == "source_url":
        kwargs["source_name"] = "Exact source"
    kwargs[argument] = value

    with pytest.raises(ValueError) as failure:
        prepare(out_path, format="h5ad", **kwargs)

    message = str(failure.value)
    assert message.startswith(f"{argument} is displayed verbatim")
    assert expected in message
    assert list(tmp_path.iterdir()) == []


def test_prepare_preserves_text_that_reads_as_it_is_stored(tmp_path):
    import anndata as ad

    out_path = tmp_path / "text.h5ad"
    prepare(
        out_path,
        format="h5ad",
        embeddings={"X_umap": _EMBEDDING.copy()},
        dataset_id="display-text",
        dataset_name="Îlot de Langerhans 细胞 🧬",
        dataset_description="",
        source_name="T  cell atlas",
        gene_expression=np.ones((3, 1)),
        var=pd.DataFrame(index=["gene-1"]),
    )
    cellucid = ad.read_h5ad(out_path).uns["cellucid"]
    assert cellucid["name"] == "Îlot de Langerhans 细胞 🧬"
    assert cellucid["description"] == ""
    assert cellucid["source"] == {"name": "T  cell atlas"}


@pytest.mark.parametrize("argument", _PROSE)
@pytest.mark.parametrize(("value", "_"), _EDGES)
def test_descriptive_text_is_written_with_white_space_at_its_ends_removed(tmp_path, argument, value, _):
    import anndata as ad

    kwargs = _prepare_kwargs()
    if argument == "citation_text":
        kwargs["citation_url"] = "https://doi.org/10.1/x"
    kwargs[argument] = value
    out_path = prepare(tmp_path / "t.h5ad", format="h5ad", **kwargs)
    cellucid = ad.read_h5ad(out_path).uns["cellucid"]
    stored = {
        "dataset_name": cellucid["name"],
        "dataset_description": cellucid.get("description"),
        "citation_text": cellucid.get("citation", {}).get("text"),
        "source_name": cellucid.get("source", {}).get("name"),
    }[argument]
    assert stored == value.strip(" \u00a0\ufeff")
