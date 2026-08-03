"""The AnnData entry points must be usable on what Scanpy actually produces.

``sc.tl.umap()`` writes ``adata.obsm['X_umap']``, which names no dimension. That
object is served as written: the key resolves to the dimension its own column
count states. Cohort metadata routinely carries a ``datetime64`` column, which
cannot be served at all, so that one has to fail where the user can act on it
and name the exact action.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.anndata_server import AnnDataServer, serve_anndata
from cellucid.jupyter import AnnDataViewer, show_anndata

N_CELLS = 6


def _scanpy_shaped_adata(*, obsm_key: str = "X_umap", dimension: int = 2) -> AnnData:
    """Build the AnnData a first Scanpy run produces."""
    rng = np.random.default_rng(17)
    return AnnData(
        X=rng.random((N_CELLS, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {
                "leiden": pd.Categorical(["0", "1"] * (N_CELLS // 2)),
                "total_counts": np.arange(N_CELLS, dtype=np.float32),
            },
            index=[f"cell-{index}" for index in range(N_CELLS)],
        ),
        var=pd.DataFrame(index=["GeneA", "GeneB"]),
        obsm={
            "X_pca": rng.random((N_CELLS, 5)).astype(np.float32),
            obsm_key: rng.random((N_CELLS, dimension)).astype(np.float32),
        },
    )


def _adapter(adata: AnnData, **kwargs: object) -> AnnDataAdapter:
    return AnnDataAdapter(
        adata,
        dataset_name="Scanpy fixture",
        dataset_id="scanpy-fixture",
        centroid_min_points=1,
        **kwargs,  # type: ignore[arg-type]
    )


@pytest.mark.parametrize("dimension", [1, 2, 3])
def test_unsuffixed_scanpy_umap_is_served_at_its_own_column_count(dimension: int) -> None:
    adata = _scanpy_shaped_adata(dimension=dimension)
    before = set(adata.obsm)

    adapter = _adapter(adata)
    try:
        assert adapter.get_dataset_identity()["embeddings"]["available_dimensions"] == [dimension]
        assert adapter.get_embedding(dimension).shape == (N_CELLS, dimension)
        # The object the caller passed is the object the caller keeps: nothing
        # is written back into obsm to make this work.
        assert set(adata.obsm) == before
        assert "X_umap_2d" not in adata.obsm
    finally:
        adapter.close()


def test_an_explicit_dimensional_key_decides_over_the_unsuffixed_one() -> None:
    adata = _scanpy_shaped_adata(dimension=3)
    adata.obsm["X_umap_2d"] = np.linspace(
        -1.0,
        1.0,
        N_CELLS * 2,
        dtype=np.float32,
    ).reshape(N_CELLS, 2)

    adapter = _adapter(adata)
    try:
        # The declared key wins outright: a writer that named dimensions meant
        # to, and the 3-column X_umap does not silently join the set.
        identity = adapter.get_dataset_identity()
        assert identity["embeddings"]["available_dimensions"] == [2]
    finally:
        adapter.close()


def test_an_unsuffixed_umap_of_an_unreadable_width_names_its_own_shape() -> None:
    adata = _scanpy_shaped_adata(dimension=2)
    adata.obsm["X_umap"] = np.zeros((N_CELLS, 10), dtype=np.float32)

    with pytest.raises(ValueError) as error:
        _adapter(adata)

    message = str(error.value)
    assert "'X_umap'" in message
    assert f"({N_CELLS}, 10)" in message
    assert "X_umap_2d" in message


def test_every_unsuffixed_embedding_candidate_is_offered_with_its_own_dimension() -> None:
    # No X_umap at all, so nothing resolves and every candidate is offered.
    adata = _scanpy_shaped_adata(obsm_key="X_tsne", dimension=2)
    adata.obsm["X_phate"] = np.zeros((N_CELLS, 3), dtype=np.float32)
    # Too wide to be an embedding, so it must not be offered as one.
    adata.obsm["X_scvi"] = np.zeros((N_CELLS, 10), dtype=np.float32)

    with pytest.raises(ValueError) as error:
        _adapter(adata)

    message = str(error.value)
    assert 'adata.obsm["X_umap_2d"] = adata.obsm[\'X_tsne\']' in message
    assert 'adata.obsm["X_umap_3d"] = adata.obsm[\'X_phate\']' in message
    assert "X_scvi" not in message.split("Fix:")[1]

    # Executing exactly the suggested statement makes the same object work.
    exec('adata.obsm["X_umap_2d"] = adata.obsm[\'X_tsne\']', {"adata": adata})  # noqa: S102
    adapter = _adapter(adata)
    try:
        assert adapter.get_embedding(2).shape == (N_CELLS, 2)
    finally:
        adapter.close()


def test_declared_vector_field_keys_are_never_offered_as_embedding_renames() -> None:
    adata = _scanpy_shaped_adata(obsm_key="velocity_umap_2d")

    with pytest.raises(ValueError) as error:
        _adapter(adata)

    message = str(error.value)
    assert "velocity_umap_2d" not in message.split("Available obsm keys:")[0]


def test_transition_drift_helper_reads_the_same_unsuffixed_key() -> None:
    from cellucid.vector_fields import add_transition_drift_to_obsm

    adata = _scanpy_shaped_adata()

    # The helper resolves X_umap exactly as the adapter does, so one object
    # cannot be servable and undriftable at the same time.
    written = add_transition_drift_to_obsm(
        adata,
        np.eye(N_CELLS, dtype=np.float32),
        normalize_rows=False,
    )
    assert written == "T_fwd_umap_2d"


def test_transition_drift_helper_gives_the_same_declaration_hint() -> None:
    from cellucid.vector_fields import add_transition_drift_to_obsm

    adata = _scanpy_shaped_adata(obsm_key="X_tsne", dimension=2)

    with pytest.raises(ValueError) as error:
        add_transition_drift_to_obsm(
            adata,
            np.eye(N_CELLS, dtype=np.float32),
            normalize_rows=False,
        )

    assert 'adata.obsm["X_umap_2d"] = adata.obsm[\'X_tsne\']' in str(error.value)


def test_unservable_obs_dtype_fails_at_construction_and_names_the_escape_hatch() -> None:
    adata = _scanpy_shaped_adata(obsm_key="X_umap_2d")
    adata.obs["collection_date"] = pd.to_datetime(["2020-01-01"] * N_CELLS)

    with pytest.raises(TypeError) as error:
        _adapter(adata)

    message = str(error.value)
    assert "'collection_date'" in message
    assert "obs_keys=" in message


def test_obs_keys_serves_the_named_columns_and_nothing_else() -> None:
    adata = _scanpy_shaped_adata(obsm_key="X_umap_2d")
    adata.obs["collection_date"] = pd.to_datetime(["2020-01-01"] * N_CELLS)

    adapter = _adapter(adata, obs_keys=["leiden", "total_counts"])
    try:
        assert adapter.get_obs_keys() == ["leiden", "total_counts"]
        manifest = adapter.get_obs_manifest()
        served = {entry[1] for entry in manifest["_categoricalFields"]}
        served |= {entry[1] for entry in manifest["_continuousFields"]}
        assert served == {"leiden", "total_counts"}
        # Only the two served indices resolve; a third payload route does not exist.
        assert adapter.get_obs_key_for_payload_index("2") is None
        with pytest.raises(KeyError, match="collection_date"):
            adapter.get_obs_categorical_codes("collection_date")
    finally:
        adapter.close()


@pytest.mark.parametrize(
    "invalid",
    ["leiden", [], ["leiden", "not_a_column"]],
)
def test_obs_keys_rejects_every_selection_that_cannot_be_served(invalid: object) -> None:
    adata = _scanpy_shaped_adata(obsm_key="X_umap_2d")
    with pytest.raises((TypeError, ValueError, KeyError)):
        _adapter(adata, obs_keys=invalid)


def test_every_public_anndata_entry_point_exposes_obs_keys() -> None:
    for entry_point in (
        AnnDataAdapter.__init__,
        AnnDataAdapter.from_file,
        AnnDataServer.__init__,
        serve_anndata,
        AnnDataViewer.__init__,
        show_anndata,
    ):
        parameters = inspect.signature(entry_point).parameters
        assert "obs_keys" in parameters, entry_point.__qualname__
        assert parameters["obs_keys"].kind is inspect.Parameter.KEYWORD_ONLY


def test_from_file_serves_only_the_named_obs_columns(tmp_path: Path) -> None:
    adata = _scanpy_shaped_adata(obsm_key="X_umap_2d")
    path = tmp_path / "scanpy.h5ad"
    adata.write_h5ad(path)

    adapter = AnnDataAdapter.from_file(
        path,
        dataset_name="Scanpy fixture",
        dataset_id="scanpy-fixture",
        centroid_min_points=1,
        obs_keys=["leiden"],
    )
    try:
        assert adapter.get_obs_keys() == ["leiden"]
    finally:
        adapter.close()


def test_cli_serve_forwards_repeated_obs_key_options() -> None:
    from cellucid.cli import create_parser

    arguments = create_parser().parse_args(
        [
            "serve",
            "data.h5ad",
            "--dataset-name",
            "Example",
            "--dataset-id",
            "example",
            "--obs-key",
            "leiden",
            "--obs-key",
            "total_counts",
        ]
    )
    assert arguments.obs_key == ["leiden", "total_counts"]
