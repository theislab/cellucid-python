"""Exact contract for export format v3: embeddings as identified items.

``dataset_identity.json`` version 3 declares any number of embeddings per
dimension under ``embeddings.items``. Each item is ordered by (dimension, id in
code-point order); its position is its payload index and the one file the
format derives for it, ``embeddings/{index}_{dimension}d.bin`` (plus ``.gz``
when the export is compressed). Every dimension the export carries has one
default item, the highest dimension is the initial one, and each item records
the 0-based source columns it was read from. Every one of those rules is
checked before the first output path is created.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cellucid
from cellucid.prepare_data import Embedding, prepare

N_CELLS = 4


def _plane(scale: float = 1.0, offset: float = 0.0) -> np.ndarray:
    return (
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]], dtype=np.float64) * scale
        + offset
    )


def _kwargs(out_dir: Path, embeddings: object, **overrides: object) -> dict[str, object]:
    kwargs: dict[str, object] = {
        "latent_space": _plane(),
        "obs": pd.DataFrame(
            {"group": pd.Categorical(["a", "a", "b", "b"])},
            index=[f"cell-{index}" for index in range(N_CELLS)],
        ),
        "embeddings": embeddings,
        "out_dir": out_dir,
        "dataset_name": "Export v3",
        "dataset_id": "export-v3",
        "created_at": "2026-01-01T00:00:00Z",
        "obs_categorical_dtype": "uint8",
        "centroid_min_points": 1,
    }
    kwargs.update(overrides)
    return kwargs


def _identity(out_dir: Path) -> dict:
    return json.loads((out_dir / "dataset_identity.json").read_text(encoding="utf-8"))


def _many_embeddings() -> dict[str, object]:
    # Supplied out of order on purpose: the export orders them itself.
    return {
        "X_umap_3d": np.column_stack([_plane(), np.arange(N_CELLS, dtype=np.float64)]),
        "X_umap": _plane(),
        "X_pca": Embedding(
            np.column_stack([_plane(4.0), _plane(2.0), np.zeros(N_CELLS)]),
            columns=[3, 1],
            label="PCA 4 and 2",
            source_key="X_pca",
        ),
        "X_line": np.linspace(0.0, 1.0, N_CELLS),
        "X_TSNE": Embedding(_plane(8.0), label="t-SNE", source_key="X_tsne"),
    }


def test_identity_is_version_3_with_the_exact_embedding_shape(tmp_path: Path) -> None:
    prepare(**_kwargs(tmp_path, _many_embeddings(), embedding_defaults={"2d": "X_umap"}))
    identity = _identity(tmp_path)

    assert identity["version"] == 3
    assert identity["cellucid_data_version"] == cellucid.__version__
    assert list(identity) == [
        "version",
        "id",
        "name",
        "description",
        "created_at",
        "cellucid_data_version",
        "stats",
        "embeddings",
        "obs_fields",
        "export_settings",
    ]
    embeddings = identity["embeddings"]
    assert list(embeddings) == ["default_dimension", "defaults", "items"]
    assert embeddings["default_dimension"] == 3
    assert embeddings["defaults"] == {"1d": "X_line", "2d": "X_umap", "3d": "X_umap_3d"}
    assert list(embeddings["defaults"]) == ["1d", "2d", "3d"]
    assert embeddings["items"] == {
        "X_line": {
            "label": "X_line",
            "dimension": 1,
            "file": "embeddings/0_1d.bin",
            "source": {"key": None, "columns": [0]},
        },
        "X_TSNE": {
            "label": "t-SNE",
            "dimension": 2,
            "file": "embeddings/1_2d.bin",
            "source": {"key": "X_tsne", "columns": [0, 1]},
        },
        "X_pca": {
            "label": "PCA 4 and 2",
            "dimension": 2,
            "file": "embeddings/2_2d.bin",
            "source": {"key": "X_pca", "columns": [3, 1]},
        },
        "X_umap": {
            "label": "X_umap",
            "dimension": 2,
            "file": "embeddings/3_2d.bin",
            "source": {"key": None, "columns": [0, 1]},
        },
        "X_umap_3d": {
            "label": "X_umap_3d",
            "dimension": 3,
            "file": "embeddings/4_3d.bin",
            "source": {"key": None, "columns": [0, 1, 2]},
        },
    }
    # Order is (dimension, code point): 'X_TSNE' sorts before 'X_pca'.
    assert list(embeddings["items"]) == ["X_line", "X_TSNE", "X_pca", "X_umap", "X_umap_3d"]
    for item in embeddings["items"].values():
        assert list(item) == ["label", "dimension", "file", "source"]
        assert list(item["source"]) == ["key", "columns"]
    assert sorted(path.name for path in (tmp_path / "embeddings").iterdir()) == [
        "0_1d.bin",
        "1_2d.bin",
        "2_2d.bin",
        "3_2d.bin",
        "4_3d.bin",
    ]


@pytest.mark.parametrize(("compression", "suffix"), [(None, ""), (1, ".gz"), (9, ".gz")])
def test_each_item_file_is_exactly_the_path_its_position_derives(
    tmp_path: Path,
    compression: int | None,
    suffix: str,
) -> None:
    prepare(
        **_kwargs(
            tmp_path,
            {"b": _plane(), "a": _plane(2.0)},
            embedding_defaults={"2d": "b"},
            compression=compression,
        )
    )
    items = _identity(tmp_path)["embeddings"]["items"]
    assert [item["file"] for item in items.values()] == [
        f"embeddings/0_2d.bin{suffix}",
        f"embeddings/1_2d.bin{suffix}",
    ]
    assert list(items) == ["a", "b"]
    assert sorted(path.name for path in (tmp_path / "embeddings").iterdir()) == [
        f"0_2d.bin{suffix}",
        f"1_2d.bin{suffix}",
    ]


def test_each_item_is_normalized_on_its_own_from_its_selected_columns(tmp_path: Path) -> None:
    wide = np.column_stack([_plane(4.0, 10.0), _plane(2.0, -3.0)])
    prepare(
        **_kwargs(
            tmp_path,
            {"X_umap": _plane(), "X_pca": Embedding(wide, columns=[3, 0])},
            embedding_defaults={"2d": "X_umap"},
        )
    )
    selected = wide[:, [3, 0]]
    mins, maxs = selected.min(axis=0), selected.max(axis=0)
    expected = ((selected - (mins + maxs) / 2) * (2.0 / float((maxs - mins).max()))).astype(
        np.float32
    )
    published = np.fromfile(tmp_path / "embeddings" / "0_2d.bin", dtype="<f4").reshape(N_CELLS, 2)
    np.testing.assert_array_equal(published, expected)
    umap = np.fromfile(tmp_path / "embeddings" / "1_2d.bin", dtype="<f4").reshape(N_CELLS, 2)
    np.testing.assert_array_equal(umap, (_plane() * 2.0 - 1.0).astype(np.float32))


def test_a_single_item_per_dimension_is_that_dimensions_default(tmp_path: Path) -> None:
    prepare(**_kwargs(tmp_path, {"X_tsne": _plane(), "X_line": np.arange(N_CELLS)}))
    embeddings = _identity(tmp_path)["embeddings"]
    assert embeddings["defaults"] == {"1d": "X_line", "2d": "X_tsne"}
    assert embeddings["default_dimension"] == 2


def test_a_dimension_with_several_items_requires_a_named_default(tmp_path: Path) -> None:
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(
        ValueError,
        match=(
            r"2 2D embeddings were supplied \(\['X_tsne', 'X_umap'\]\), so "
            r"embedding_defaults must name the one shown first"
        ),
    ):
        prepare(**_kwargs(out_dir, {"X_umap": _plane(), "X_tsne": _plane(2.0)}))
    assert not out_dir.exists()


@pytest.mark.parametrize(
    ("defaults", "error", "message"),
    [
        ({"2d": "X_line"}, ValueError, r"'X_line', which is not a 2D embedding"),
        ({"2d": "missing"}, ValueError, r"'missing', which is not a 2D embedding"),
        ({"3d": "X_umap"}, ValueError, r"names a 3d default, and no 3D embedding was supplied"),
        ({"4d": "X_umap"}, ValueError, r"key '4d' must be '1d', '2d', or '3d'"),
        ({"2D": "X_umap"}, ValueError, r"key '2D' must be '1d', '2d', or '3d'"),
        ({"2d": 7}, ValueError, r"is 7, which is not a 2D embedding"),
        ("X_umap", TypeError, r"embedding_defaults must map '1d', '2d', or '3d'"),
        ([("2d", "X_umap")], TypeError, r"embedding_defaults must map"),
    ],
)
def test_invalid_embedding_defaults_are_rejected_before_output_mutation(
    tmp_path: Path,
    defaults: object,
    error: type[Exception],
    message: str,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(error, match=message):
        prepare(
            **_kwargs(
                out_dir,
                {"X_umap": _plane(), "X_tsne": _plane(2.0), "X_line": np.arange(N_CELLS)},
                embedding_defaults=defaults,
            )
        )
    assert not out_dir.exists()


def test_a_named_default_may_confirm_a_single_item(tmp_path: Path) -> None:
    prepare(**_kwargs(tmp_path, {"X_umap": _plane()}, embedding_defaults={"2d": "X_umap"}))
    assert _identity(tmp_path)["embeddings"]["defaults"] == {"2d": "X_umap"}


@pytest.mark.parametrize(
    ("embeddings", "error", "message"),
    [
        (None, ValueError, r"embeddings is required"),
        ({}, ValueError, r"embeddings must contain at least one embedding"),
        ("X_umap", TypeError, r"embeddings must be a mapping"),
        ([_plane()], TypeError, r"embeddings must be a mapping"),
        ({"": _plane()}, ValueError, r"Embedding identifier at position 0 must be non-empty"),
        ({7: _plane()}, TypeError, r"Embedding identifier at position 0 must be a string"),
        ({" X_umap": _plane()}, ValueError, r"displayed verbatim"),
        ({"X_umap": None}, ValueError, r"Embedding 'X_umap' has no values"),
        ({"X_umap": np.zeros((N_CELLS, 2, 2))}, ValueError, r"2D array or a 1D vector"),
        ({"X_umap": np.array([["a", "b"]] * N_CELLS)}, TypeError, r"real numeric"),
        ({"X_umap": np.ones((N_CELLS, 2), dtype=bool)}, TypeError, r"real numeric"),
        (
            {"X_umap": _plane(), "X_tsne": _plane()[:3]},
            ValueError,
            r"same number of cells: 'X_umap' has 4 and 'X_tsne' has 3",
        ),
        ({"X_umap": np.ones((N_CELLS, 2))}, ValueError, r"no coordinate variation"),
        (
            {"X_pca": np.zeros((N_CELLS, 5))},
            ValueError,
            r"'X_pca' has 5 columns, and an embedding plots 1, 2, or 3 of them\. "
            r"Pass Embedding\(values, columns=\[\.\.\.\]\)",
        ),
        (
            {"X_umap": Embedding(_plane(), label=" UMAP")},
            ValueError,
            r"Embedding 'X_umap' label is displayed verbatim",
        ),
        ({"X_umap": Embedding(_plane(), label="")}, ValueError, r"label must be a non-empty"),
        (
            {"X_umap": Embedding(_plane(), source_key="")},
            ValueError,
            r"source_key must be a non-empty string",
        ),
    ],
)
def test_invalid_embeddings_are_rejected_before_output_mutation(
    tmp_path: Path,
    embeddings: object,
    error: type[Exception],
    message: str,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(error, match=message):
        prepare(**_kwargs(out_dir, embeddings))
    assert not out_dir.exists()


@pytest.mark.parametrize(
    ("columns", "error", "message"),
    [
        ([1, 1], ValueError, r"'X_pca': Each axis needs a different column\."),
        ([0, 5], ValueError, r"'X_pca': Embedding column 5 does not exist; the array has columns 0 to 4"),
        ([-1, 0], ValueError, r"column -1 does not exist"),
        ([0, 1.0], ValueError, r"Embedding column 1\.0 must be an integer"),
        ([0, True], ValueError, r"Embedding column True must be an integer"),
        ([], ValueError, r"columns must name 1, 2, or 3 columns, one per axis; got 0"),
        ([0, 1, 2, 3], ValueError, r"one per axis; got 4"),
        ("01", TypeError, r"columns must be a sequence of 0-based integers"),
        (3, TypeError, r"columns must be a sequence of 0-based integers"),
    ],
)
def test_invalid_embedding_columns_are_rejected_before_output_mutation(
    tmp_path: Path,
    columns: object,
    error: type[Exception],
    message: str,
) -> None:
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(error, match=message):
        prepare(**_kwargs(out_dir, {"X_pca": Embedding(np.zeros((N_CELLS, 5)), columns=columns)}))  # type: ignore[arg-type]
    assert not out_dir.exists()


def test_columns_wider_than_the_array_are_too_narrow(tmp_path: Path) -> None:
    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(
        ValueError,
        match=r"'X_line': A 3D embedding needs at least 3 columns; this array has 1",
    ):
        prepare(**_kwargs(out_dir, {"X_line": Embedding(np.arange(N_CELLS), columns=[0, 1, 2])}))
    assert not out_dir.exists()


def test_only_the_selected_columns_must_be_finite(tmp_path: Path) -> None:
    wide = np.column_stack([_plane(), np.full(N_CELLS, np.nan), _plane(3.0)])
    prepare(**_kwargs(tmp_path, {"X_pca": Embedding(wide, columns=[0, 4])}))
    assert _identity(tmp_path)["embeddings"]["items"]["X_pca"]["source"]["columns"] == [0, 4]

    out_dir = tmp_path / "must-not-exist"
    with pytest.raises(ValueError, match=r"Embedding 'X_pca' must contain only finite values"):
        prepare(**_kwargs(out_dir, {"X_pca": Embedding(wide, columns=[2, 0])}))
    assert not out_dir.exists()


def test_columns_may_reorder_and_subset_a_narrow_array(tmp_path: Path) -> None:
    source = np.column_stack([_plane(), np.arange(N_CELLS, dtype=np.float64)])
    prepare(**_kwargs(tmp_path, {"X_umap_3d": Embedding(source, columns=[2, 0])}))
    item = _identity(tmp_path)["embeddings"]["items"]["X_umap_3d"]
    assert item["dimension"] == 2
    assert item["file"] == "embeddings/0_2d.bin"
    assert item["source"]["columns"] == [2, 0]


def test_categorical_centroids_are_keyed_by_every_item_in_item_order(tmp_path: Path) -> None:
    prepare(**_kwargs(tmp_path, _many_embeddings(), embedding_defaults={"2d": "X_umap"}))
    manifest = json.loads((tmp_path / "obs_manifest.json").read_text(encoding="utf-8"))
    identity = _identity(tmp_path)

    assert manifest["_format"] == "compact_v2"
    (field,) = manifest["_categoricalFields"]
    centroids_by_embedding = field[5]
    assert list(centroids_by_embedding) == list(identity["embeddings"]["items"])
    for embedding_id, centroids in centroids_by_embedding.items():
        dimension = identity["embeddings"]["items"][embedding_id]["dimension"]
        assert [centroid["category"] for centroid in centroids] == ["a", "b"]
        for centroid in centroids:
            assert list(centroid) == ["category", "position", "n_points"]
            assert len(centroid["position"]) == dimension


def test_centroids_are_measured_on_each_items_own_coordinates(tmp_path: Path) -> None:
    prepare(
        **_kwargs(
            tmp_path,
            {"X_umap": _plane(), "X_flipped": _plane()[:, ::-1] * 3.0},
            embedding_defaults={"2d": "X_umap"},
        )
    )
    manifest = json.loads((tmp_path / "obs_manifest.json").read_text(encoding="utf-8"))
    centroids = manifest["_categoricalFields"][0][5]
    # Group 'a' holds cells 0 and 1: (0,0),(1,0) normalized to (-1,-1),(1,-1).
    assert centroids["X_umap"][0]["position"] == [0.0, -1.0]
    # In the flipped layout the same cells sit at (-1,-1),(-1,1).
    assert centroids["X_flipped"][0]["position"] == [-1.0, 0.0]


def test_none_centroid_quantile_keeps_one_empty_list_per_item(tmp_path: Path) -> None:
    prepare(
        **_kwargs(
            tmp_path,
            {"X_umap": _plane(), "X_tsne": _plane(2.0)},
            embedding_defaults={"2d": "X_umap"},
            centroid_outlier_quantile=None,
        )
    )
    manifest = json.loads((tmp_path / "obs_manifest.json").read_text(encoding="utf-8"))
    assert manifest["_categoricalFields"][0][5] == {"X_tsne": [], "X_umap": []}


def test_the_embedding_type_is_exported_from_the_package() -> None:
    assert cellucid.Embedding is Embedding
    assert "Embedding" in cellucid.__all__
    assert Embedding(_plane()).columns is None
    with pytest.raises(AttributeError):
        Embedding(_plane()).label = "frozen"  # type: ignore[misc]
