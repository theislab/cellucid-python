"""Every user-supplied string the viewer prints must read as the value it stores.

A dataset name, a description, a citation and a source line are drawn verbatim
in the viewer. A character with no glyph therefore changes the stored value
without changing the picture, so ``prepare()`` rejects such text instead of
trimming it: trimming rewrites text nobody asked to change.
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


@pytest.mark.parametrize(
    "argument",
    [
        "dataset_name",
        "dataset_description",
        "citation_text",
        "citation_url",
        "source_name",
        "source_url",
    ],
)
@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (" padded ", "starts with U+0020 SPACE"),
        ("padded ", "ends with U+00A0 NO-BREAK SPACE"),
        ("hid​den", "contains U+200B ZERO WIDTH SPACE"),
        ("cont\x01rol", "contains U+0001 (control character)"),
        ("Li\ud800ver", "contains U+D800 (unpaired surrogate)"),
        ("﻿Liver", "contains U+FEFF ZERO WIDTH NO-BREAK SPACE"),
    ],
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
