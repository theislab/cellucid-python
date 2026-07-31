"""Every user-supplied string the viewer prints must read as the value it stores.

A category label, a dataset name, a description, and a source line are all drawn
verbatim -- in the legend, the field selector, and every exported figure. A
character with no glyph therefore changes the stored value without changing the
picture, which is how a published export ended up with the eight organs
``'Bone_marrow '`` through ``'Yolk_sac '`` beside the one clean
``'Mesenteric_lymph_node'``.

The exporter rejects those strings instead of trimming them: trimming rewrites an
annotation nobody asked to change, and merges ``'Liver '`` into a separate
``'Liver'`` category, moving cells between them without saying so.
"""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.prepare_data import generate_datasets_manifest, prepare

_EMBEDDING = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float32)

# One representative of every rejected class, with the exact phrase the message
# must use to name it.
UNDISPLAYABLE_LABELS = [
    pytest.param(["Liver ", "Gut ", "Liver "], "ends with U+0020 SPACE", id="trailing-space"),
    pytest.param([" Liver", "Gut", "Gut"], "starts with U+0020 SPACE", id="leading-space"),
    pytest.param(["   ", "Gut", "Gut"], "starts with U+0020 SPACE", id="whitespace-only"),
    pytest.param(["", "Gut", "Gut"], "'' is empty", id="empty"),
    pytest.param(
        ["Liver ", "Gut", "Gut"],
        "ends with U+00A0 NO-BREAK SPACE",
        id="nbsp",
    ),
    pytest.param(
        ["Liver　", "Gut", "Gut"],
        "ends with U+3000 IDEOGRAPHIC SPACE",
        id="ideographic-space",
    ),
    pytest.param(
        ["Li\tver", "Gut", "Gut"],
        "contains U+0009 (control character)",
        id="tab",
    ),
    pytest.param(
        ["Li\nver", "Gut", "Gut"],
        "contains U+000A (control character)",
        id="newline",
    ),
    pytest.param(
        ["Liver\x85", "Gut", "Gut"],
        "contains U+0085 (control character)",
        id="c1-control",
    ),
    pytest.param(
        ["Liver​", "Gut", "Gut"],
        "contains U+200B ZERO WIDTH SPACE",
        id="zero-width-space",
    ),
    pytest.param(
        ["﻿Liver", "Gut", "Gut"],
        "contains U+FEFF ZERO WIDTH NO-BREAK SPACE",
        id="byte-order-mark",
    ),
    pytest.param(
        ["Liver⁠", "Gut", "Gut"],
        "contains U+2060 WORD JOINER",
        id="word-joiner",
    ),
]


def _prepare_kwargs(out_dir, obs):
    return {
        "latent_space": _EMBEDDING.copy(),
        "obs": obs,
        "X_umap_2d": _EMBEDDING.copy(),
        "out_dir": out_dir,
        "dataset_id": "display-text",
        "dataset_name": "Display text",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
    }


def _adata(obs: pd.DataFrame) -> ad.AnnData:
    adata = ad.AnnData(
        X=np.array([[0.0], [1.0], [2.0]], dtype=np.float32),
        obs=obs.set_axis(pd.Index(["cell-1", "cell-2", "cell-3"], dtype=object)),
        var=pd.DataFrame(index=pd.Index(["Gene_A"], dtype=object)),
    )
    adata.obsm["X_umap_2d"] = _EMBEDDING.copy()
    adata.obsm["X_latent"] = _EMBEDDING.copy()
    return adata


@pytest.mark.parametrize(("labels", "expected"), UNDISPLAYABLE_LABELS)
def test_prepare_rejects_undisplayable_category_labels_before_output_mutation(
    tmp_path,
    labels,
    expected,
):
    out_dir = tmp_path / "must-not-exist"
    obs = pd.DataFrame({"organ": labels})

    with pytest.raises(ValueError) as failure:
        prepare(**_prepare_kwargs(out_dir, obs))

    message = str(failure.value)
    assert "Categorical field 'organ'" in message
    assert expected in message
    assert not out_dir.exists()


@pytest.mark.parametrize(("labels", "expected"), UNDISPLAYABLE_LABELS)
def test_adapter_rejects_undisplayable_category_labels(labels, expected):
    adapter = AnnDataAdapter(
        _adata(pd.DataFrame({"organ": labels})),
        latent_key="X_latent",
        dataset_name="Display text",
        dataset_id="display-text",
    )

    with pytest.raises(ValueError) as failure:
        adapter.get_obs_manifest()

    message = str(failure.value)
    assert "Categorical field 'organ'" in message
    assert expected in message


def test_prepare_rejects_the_published_suo_organ_labels(tmp_path):
    """The exact labels the published ``suo`` export shipped."""
    out_dir = tmp_path / "must-not-exist"
    published = [
        "Bone_marrow ",
        "Gut ",
        "Kidney ",
        "Liver ",
        "Mesenteric_lymph_node",
        "Skin ",
        "Spleen ",
        "Thymus ",
        "Yolk_sac ",
    ]
    obs = pd.DataFrame({"organ": pd.Categorical(published * 3, categories=published)})
    kwargs = _prepare_kwargs(out_dir, obs)
    kwargs["latent_space"] = np.tile(_EMBEDDING, (9, 1))
    kwargs["X_umap_2d"] = np.tile(_EMBEDDING, (9, 1))

    with pytest.raises(ValueError) as failure:
        prepare(**kwargs)

    message = str(failure.value)
    # Every offending label is named at once; the one clean label is not.
    assert "8 labels" in message
    for label in published[:4]:
        assert repr(label) in message
    assert "'Mesenteric_lymph_node'" not in message
    assert not out_dir.exists()


def test_prepare_rejects_a_padded_label_colliding_with_its_clean_twin(tmp_path):
    """'Liver' and 'Liver ' are distinct categories that look identical."""
    out_dir = tmp_path / "must-not-exist"
    obs = pd.DataFrame({"organ": ["Liver", "Liver ", "Liver"]})

    with pytest.raises(ValueError) as failure:
        prepare(**_prepare_kwargs(out_dir, obs))

    message = str(failure.value)
    assert "'Liver '" in message
    assert "ends with U+0020 SPACE" in message
    # The remedy must not read as an instruction to merge the two categories.
    assert "can merge two categories into one" in message
    assert not out_dir.exists()


def test_prepare_rejects_labels_a_whitespace_collapsing_renderer_draws_alike(tmp_path):
    out_dir = tmp_path / "must-not-exist"
    obs = pd.DataFrame({"organ": ["T cell", "T  cell", "T cell"]})

    with pytest.raises(ValueError) as failure:
        prepare(**_prepare_kwargs(out_dir, obs))

    message = str(failure.value)
    assert "'T  cell'" in message
    assert "'T cell'" in message
    assert "drawn identically" in message
    assert not out_dir.exists()


def test_prepare_names_the_field_and_the_exact_repair(tmp_path):
    out_dir = tmp_path / "must-not-exist"
    obs = pd.DataFrame({"organ": ["Liver ", "Gut", "Gut"]})

    with pytest.raises(ValueError) as failure:
        prepare(**_prepare_kwargs(out_dir, obs))

    message = str(failure.value)
    assert "the legend, the field selector, and exported figures" in message
    assert "Cellucid does not clean them for you" in message
    assert "obs['organ'] = obs['organ'].astype('string').str.strip()" in message
    assert "obs['organ'] = obs['organ'].astype('category')" in message


def test_prepare_preserves_every_label_that_reads_as_it_is_stored(tmp_path):
    out_dir = tmp_path / "generation"
    obs = pd.DataFrame(
        {
            "organ": ["Îlot de Langerhans", "T cell", "细胞 🧬"],
            "case": ["Liver", "liver", "LIVER"],
            "flag": [True, False, True],
            "cluster": pd.Categorical([1, 2, 1]),
        }
    )

    prepare(**_prepare_kwargs(out_dir, obs))

    manifest = json.loads((out_dir / "obs_manifest.json").read_text(encoding="utf-8"))
    categories = {field[0]: field[1] for field in manifest["_categoricalFields"]}
    assert categories["organ"] == ["T cell", "Îlot de Langerhans", "细胞 🧬"]
    # Case-only differences are visible on screen, so they stay two categories.
    assert categories["case"] == ["LIVER", "Liver", "liver"]
    assert categories["flag"] == [False, True]
    assert categories["cluster"] == [1, 2]


@pytest.mark.parametrize(
    "argument",
    ["dataset_name", "dataset_description", "source_name", "source_url", "source_citation"],
)
@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (" padded ", "starts with U+0020 SPACE"),
        ("padded ", "ends with U+00A0 NO-BREAK SPACE"),
        ("hid​den", "contains U+200B ZERO WIDTH SPACE"),
        ("cont\x01rol", "contains U+0001 (control character)"),
    ],
)
def test_prepare_rejects_undisplayable_identity_text_before_output_mutation(
    tmp_path,
    argument,
    value,
    expected,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_kwargs(out_dir, pd.DataFrame(index=["a", "b", "c"]))
    if argument in {"source_url", "source_citation"}:
        kwargs["source_name"] = "Exact source"
    kwargs[argument] = value

    with pytest.raises(ValueError) as failure:
        prepare(**kwargs)

    message = str(failure.value)
    assert message.startswith(f"{argument} is displayed verbatim")
    assert expected in message
    assert not out_dir.exists()
    assert list(tmp_path.iterdir()) == []


def test_datasets_manifest_rejects_an_undisplayable_published_description(tmp_path):
    exports = tmp_path / "exports"
    dataset = exports / "study"
    dataset.mkdir(parents=True)
    (dataset / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": "study",
                "name": "Study",
                "description": "Padded description ",
                "stats": {"n_cells": 3, "n_genes": 1},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError) as failure:
        generate_datasets_manifest(exports_dir=exports, default_dataset="study")

    assert "description is displayed verbatim" in str(failure.value)
    assert "ends with U+0020 SPACE" in str(failure.value)
    assert not (exports / "datasets.json").exists()
