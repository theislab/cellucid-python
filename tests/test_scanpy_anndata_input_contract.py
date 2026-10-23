"""The AnnData entry points must be usable on what Scanpy actually produces.

``sc.tl.umap()`` writes ``adata.obsm['X_umap']``, which names no dimension. That
object is served as written: the key resolves to the dimension its own column
count states. An object with no UMAP at all is still served: Cellucid makes a
deterministic best guess from the other ``obsm`` arrays and says so. Cohort metadata routinely carries a ``datetime64`` column, which
cannot be served at all, so that one has to fail where the user can act on it
and name the exact action.
"""

from __future__ import annotations

import inspect

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

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


def test_transition_drift_helper_reads_the_same_unsuffixed_key() -> None:
    from cellucid.vector_fields import add_transition_drift_to_obsm

    adata = _scanpy_shaped_adata()

    # The helper reads X_umap, the key sc.tl.umap() writes, and names the
    # drift for that basis, the key the viewer binds to X_umap.
    written = add_transition_drift_to_obsm(
        adata,
        np.eye(N_CELLS, dtype=np.float32),
        normalize_rows=False,
    )
    assert written == "T_fwd_umap"


def test_transition_drift_helper_names_the_umap_keys_it_reads() -> None:
    from cellucid.vector_fields import add_transition_drift_to_obsm

    adata = _scanpy_shaped_adata(obsm_key="X_tsne", dimension=2)

    # A UMAP drift is computed in UMAP coordinates, so another embedding is
    # not substituted: the helper names what it reads and what obsm holds.
    with pytest.raises(ValueError) as error:
        add_transition_drift_to_obsm(
            adata,
            np.eye(N_CELLS, dtype=np.float32),
            normalize_rows=False,
        )

    message = str(error.value)
    assert "obsm['X_umap']" in message
    assert "'X_tsne'" in message


def test_every_public_anndata_entry_point_exposes_obs_keys() -> None:
    for entry_point in (
        AnnDataServer.__init__,
        serve_anndata,
        AnnDataViewer.__init__,
        show_anndata,
    ):
        parameters = inspect.signature(entry_point).parameters
        assert "obs_keys" in parameters, entry_point.__qualname__
        assert parameters["obs_keys"].kind is inspect.Parameter.KEYWORD_ONLY


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
