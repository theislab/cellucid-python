import json

import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data import generate_datasets_manifest, prepare


def _prepare_identity_kwargs(out_dir):
    embedding = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ],
        dtype=np.float32,
    )
    return {
        "latent_space": embedding.copy(),
        "obs": pd.DataFrame(index=["cell-1", "cell-2", "cell-3"]),
        "X_umap_2d": embedding,
        "out_dir": out_dir,
        "dataset_id": "identity-contract",
        "dataset_name": "Identity contract",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
    }


def _prepare_obs_export(tmp_path, obs, *, force):
    prepare(
        latent_space=np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
            ],
            dtype=np.float32,
        ),
        obs=obs,
        X_umap_2d=np.array(
            [
                [0.0, 0.0],
                [1.0, 0.0],
                [0.0, 1.0],
            ],
            dtype=np.float32,
        ),
        out_dir=tmp_path,
        dataset_id="obs-contract",
        dataset_name="Observation contract",
        obs_categorical_dtype="uint16",
        centroid_min_points=1,
        force=force,
    )


def test_identity_obs_fields_follow_emitted_compact_manifest_order(tmp_path):
    obs = pd.DataFrame(
        {
            "group": pd.Categorical(["A", "A", "B"]),
            "score": [0.1, 0.2, 0.3],
            "selected": [True, False, True],
            "quality": [1.0, 2.0, 3.0],
        }
    )

    _prepare_obs_export(tmp_path, obs, force=True)

    manifest = json.loads((tmp_path / "obs_manifest.json").read_text(encoding="utf-8"))
    identity = json.loads((tmp_path / "dataset_identity.json").read_text(encoding="utf-8"))

    manifest_order = [
        *({"key": field[0], "kind": "continuous"} for field in manifest["_continuousFields"]),
        *(
            {
                "key": field[0],
                "kind": "category",
                "n_categories": len(field[1]),
            }
            for field in manifest["_categoricalFields"]
        ),
    ]
    assert identity["obs_fields"] == manifest_order
    assert identity["obs_fields"] == [
        {"key": "score", "kind": "continuous"},
        {"key": "quality", "kind": "continuous"},
        {"key": "group", "kind": "category", "n_categories": 2},
        {"key": "selected", "kind": "category", "n_categories": 2},
    ]
    assert identity["stats"]["n_obs_fields"] == 4
    assert identity["stats"]["n_continuous_fields"] == 2
    assert identity["stats"]["n_categorical_fields"] == 2


def test_existing_obs_manifest_must_match_requested_fields(tmp_path):
    original_obs = pd.DataFrame(
        {
            "group": pd.Categorical(["A", "A", "B"]),
            "score": [0.1, 0.2, 0.3],
        }
    )
    changed_obs = pd.DataFrame(
        {
            "group": pd.Categorical(["A", "A", "B"]),
            "different_score": [0.1, 0.2, 0.3],
        }
    )
    _prepare_obs_export(tmp_path, original_obs, force=True)

    with pytest.raises(FileExistsError, match="non-empty export directory"):
        _prepare_obs_export(tmp_path, changed_obs, force=False)


@pytest.mark.parametrize("omitted_key", ["dataset_name", "dataset_id"])
def test_prepare_requires_explicit_dataset_identity_before_output_mutation(
    tmp_path,
    omitted_key,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    del kwargs[omitted_key]

    with pytest.raises(TypeError, match=omitted_key):
        prepare(**kwargs)

    assert not out_dir.exists()


@pytest.mark.parametrize(
    "dataset_id",
    [
        "",
        "   ",
        " padded ",
        "data/set",
        "data\\set",
        "CON",
        "trailing.",
        "ümlaut",
        "a" * 181,
    ],
)
def test_prepare_rejects_nonportable_dataset_id_before_output_mutation(
    tmp_path,
    dataset_id,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["dataset_id"] = dataset_id

    with pytest.raises((TypeError, ValueError), match="dataset_id"):
        prepare(**kwargs)

    assert not out_dir.exists()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    "dataset_name",
    [
        "",
        "   ",
        " padded",
        "padded ",
        "\tTabbed",
        "Control\x7f",
        "Control\x80",
        "Control\x85",
        "Control\x9f",
    ],
)
def test_prepare_rejects_inexact_dataset_name_before_output_mutation(
    tmp_path,
    dataset_name,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["dataset_name"] = dataset_name

    with pytest.raises((TypeError, ValueError), match="dataset_name"):
        prepare(**kwargs)

    assert not out_dir.exists()
    assert list(tmp_path.iterdir()) == []


def test_prepare_preserves_unicode_dataset_name_with_portable_id(tmp_path):
    out_dir = tmp_path / "generation"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["dataset_name"] = "细胞 atlas 🧬"

    prepare(**kwargs)

    identity = json.loads((out_dir / "dataset_identity.json").read_text(encoding="utf-8"))
    assert identity["id"] == "identity-contract"
    assert identity["name"] == "细胞 atlas 🧬"


@pytest.mark.parametrize("invalid_description", [7, True, [], {}])
def test_prepare_rejects_non_string_description_before_output_mutation(
    tmp_path,
    invalid_description,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["dataset_description"] = invalid_description

    with pytest.raises(TypeError, match="dataset_description"):
        prepare(**kwargs)

    assert not out_dir.exists()
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize(
    ("source_overrides", "message"),
    [
        ({"source_url": "https://example.test"}, "source_name"),
        ({"source_citation": "Reference"}, "source_name"),
        ({"source_name": 7}, "source_name"),
        ({"source_name": ""}, "source_name"),
        ({"source_name": "   "}, "source_name"),
        ({"source_name": "Source", "source_url": 7}, "source_url"),
        ({"source_name": "Source", "source_url": ""}, "source_url"),
        ({"source_name": "Source", "source_citation": 7}, "source_citation"),
        ({"source_name": "Source", "source_citation": ""}, "source_citation"),
    ],
)
def test_prepare_rejects_partial_or_non_string_source_before_output_mutation(
    tmp_path,
    source_overrides,
    message,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs.update(source_overrides)

    with pytest.raises((TypeError, ValueError), match=message):
        prepare(**kwargs)

    assert not out_dir.exists()
    assert list(tmp_path.iterdir()) == []


def test_prepare_emits_one_exact_optional_source_object(tmp_path):
    out_dir = tmp_path / "generation"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs.update(
        {
            "dataset_description": "",
            "source_name": "Exact source",
            "source_url": "https://example.test/study",
            "source_citation": "Exact citation",
        }
    )
    prepare(**kwargs)

    identity = json.loads((out_dir / "dataset_identity.json").read_text(encoding="utf-8"))
    assert identity["description"] == ""
    assert identity["source"] == {
        "name": "Exact source",
        "url": "https://example.test/study",
        "citation": "Exact citation",
    }


@pytest.mark.parametrize("invalid_value", [np.nan, np.inf, -np.inf])
def test_prepare_rejects_nonfinite_embedding_before_output_mutation(
    tmp_path,
    invalid_value,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    embedding = kwargs["X_umap_2d"].copy()
    embedding[1, 0] = invalid_value
    kwargs["X_umap_2d"] = embedding

    with pytest.raises(ValueError, match=r"X_umap_2d.*finite"):
        prepare(**kwargs)

    assert not out_dir.exists()


def test_prepare_rejects_degenerate_embedding_before_output_mutation(tmp_path):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["X_umap_2d"] = np.ones((3, 2), dtype=np.float32)

    with pytest.raises(ValueError, match=r"X_umap_2d.*no coordinate variation"):
        prepare(**kwargs)

    assert not out_dir.exists()


@pytest.mark.parametrize(
    "invalid_value",
    [np.nan, np.inf, -np.inf, np.finfo(np.float64).max],
)
def test_prepare_rejects_invalid_continuous_obs_before_output_mutation(
    tmp_path,
    invalid_value,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["obs"] = pd.DataFrame(
        {
            "score": np.array(
                [0.25, invalid_value, 0.75],
                dtype=np.float64,
            )
        }
    )

    with pytest.raises(ValueError, match=r"score.*finite|score.*float32"):
        prepare(**kwargs)

    assert not out_dir.exists()


@pytest.mark.parametrize(
    "invalid_key",
    ["bad key", "bad/key", "_leading", "trailing.", "CON", "a" * 181],
)
def test_prepare_rejects_nonportable_obs_keys_before_output_mutation(
    tmp_path,
    invalid_key,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["obs"] = pd.DataFrame({invalid_key: [1.0, 2.0, 3.0]})

    with pytest.raises(ValueError, match="ASCII|reserved Windows|end with"):
        prepare(**kwargs)

    assert not out_dir.exists()


def test_prepare_rejects_case_colliding_obs_keys_before_output_mutation(tmp_path):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["obs"] = pd.DataFrame(
        {
            "Field": [1.0, 2.0, 3.0],
            "field": [4.0, 5.0, 6.0],
        }
    )

    with pytest.raises(ValueError, match="case-insensitively"):
        prepare(**kwargs)

    assert not out_dir.exists()


@pytest.mark.parametrize(
    "gene_ids",
    [
        ["Gene_A", "bad gene"],
        ["Gene", "gene"],
        ["Gene_A", "NUL.txt"],
    ],
)
def test_prepare_rejects_nonportable_or_colliding_gene_ids_before_mutation(
    tmp_path,
    gene_ids,
):
    out_dir = tmp_path / "must-not-exist"
    kwargs = _prepare_identity_kwargs(out_dir)
    kwargs["var"] = pd.DataFrame(index=gene_ids)
    kwargs["gene_expression"] = np.array(
        [
            [0.0, 1.0],
            [2.0, 3.0],
            [4.0, 5.0],
        ],
        dtype=np.float32,
    )

    with pytest.raises(ValueError, match="ASCII|reserved Windows|case-insensitively"):
        prepare(**kwargs)

    assert not out_dir.exists()


@pytest.mark.parametrize(
    ("identity", "invalid_field"),
    [
        ({"version": 2, "name": "Declared name"}, "dataset_id"),
        ({"version": 2, "id": "declared-id"}, "dataset_name"),
        ({"version": 2, "id": 7, "name": "Declared name"}, "dataset_id"),
        ({"version": 2, "id": "data/set", "name": "Declared name"}, "dataset_id"),
        ({"version": 2, "id": " padded ", "name": "Declared name"}, "dataset_id"),
        ({"version": 2, "id": "declared-id", "name": ""}, "dataset_name"),
        ({"version": 2, "id": "declared-id", "name": "   "}, "dataset_name"),
        ({"version": 2, "id": "declared-id", "name": " padded"}, "dataset_name"),
    ],
)
def test_dataset_manifest_rejects_malformed_identity_without_filename_fabrication(
    tmp_path,
    identity,
    invalid_field,
):
    dataset_dir = tmp_path / "must-not-become-identity"
    dataset_dir.mkdir()
    (dataset_dir / "dataset_identity.json").write_text(
        json.dumps(identity),
        encoding="utf-8",
    )

    with pytest.raises((TypeError, ValueError), match=invalid_field):
        generate_datasets_manifest(tmp_path, default_dataset="declared-id")

    assert not (tmp_path / "datasets.json").exists()


@pytest.mark.parametrize(
    "default_dataset",
    ["data/set", " padded ", "CON", "trailing.", "ümlaut"],
)
def test_dataset_manifest_rejects_nonportable_default_without_mutation(
    tmp_path,
    default_dataset,
):
    dataset_dir = tmp_path / "declared-id"
    dataset_dir.mkdir()
    (dataset_dir / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": "declared-id",
                "name": "Declared name",
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="default_dataset"):
        generate_datasets_manifest(
            tmp_path,
            default_dataset=default_dataset,
        )

    assert not (tmp_path / "datasets.json").exists()
