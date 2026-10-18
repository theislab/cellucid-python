"""The ``cellucid serve`` flags that choose what an AnnData source plots and reads.

``--embedding-1d/2d/3d KEY[:c0,c1]`` names a slot's obsm array and optionally
its 0-based columns, ``--knn KEY`` the obsp neighbour graph, ``--gene-id-column``
the var column that names genes, and ``--expression X|raw|layer:NAME`` the
matrix values are read from. None of them means anything for a prepared export,
which fixed all of it when it was written, so each is refused there by name.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest import mock

import pytest

from cellucid.cli import _parse_embedding_flag, _run_serve, create_parser


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("X_tsne", "X_tsne"),
        ("X_pca:0,1", ("X_pca", [0, 1])),
        ("X_pca:4,0,2", ("X_pca", [4, 0, 2])),
        ("X_pca:7", ("X_pca", [7])),
        # A key may contain ':'; only a trailing integer list is read as columns.
        ("scvi:latent", "scvi:latent"),
        ("tool:run:3,1", ("tool:run", [3, 1])),
        ("X_pca:", "X_pca:"),
        (":0,1", ":0,1"),
        ("X_pca:0,,1", "X_pca:0,,1"),
        ("X_pca:-1,0", "X_pca:-1,0"),
    ],
)
def test_embedding_flag_grammar(value: str, expected: object) -> None:
    assert _parse_embedding_flag(value, flag="--embedding-2d") == expected


def test_an_empty_embedding_flag_names_the_flag() -> None:
    with pytest.raises(ValueError, match="--embedding-3d needs an obsm key"):
        _parse_embedding_flag("", flag="--embedding-3d")


def _zarr_store(tmp_path: Path) -> Path:
    store = tmp_path / "fixture.zarr"
    store.mkdir()
    (store / ".zgroup").write_text('{"zarr_format":2}', encoding="utf-8")
    (store / ".zattrs").write_text("{}", encoding="utf-8")
    return store


def test_every_mapping_flag_reaches_serve_anndata(tmp_path: Path) -> None:
    store = _zarr_store(tmp_path)
    args = create_parser().parse_args(
        [
            "serve",
            str(store),
            "--quiet",
            "--dataset-name",
            "Mapped",
            "--dataset-id",
            "mapped",
            "--embedding-1d",
            "X_line",
            "--embedding-2d",
            "X_pca:3,1",
            "--embedding-3d",
            "X_umap_3d",
            "--knn",
            "connectivities",
            "--gene-id-column",
            "gene_symbols",
            "--expression",
            "layer:log1p",
        ]
    )
    with mock.patch("cellucid.anndata_server.serve_anndata") as serve:
        _run_serve(args)
    kwargs = serve.call_args.kwargs
    assert kwargs["embedding_keys"] == {
        "1d": "X_line",
        "2d": ("X_pca", [3, 1]),
        "3d": "X_umap_3d",
    }
    assert kwargs["knn_key"] == "connectivities"
    assert kwargs["gene_id_column"] == "gene_symbols"
    assert kwargs["expression_source"] == "layer:log1p"
    assert "serve_connectivity" not in kwargs


def test_omitted_mapping_flags_select_the_defaults(tmp_path: Path) -> None:
    store = _zarr_store(tmp_path)
    args = create_parser().parse_args(
        ["serve", str(store), "--quiet", "--dataset-name", "D", "--dataset-id", "d"]
    )
    with mock.patch("cellucid.anndata_server.serve_anndata") as serve:
        _run_serve(args)
    kwargs = serve.call_args.kwargs
    assert kwargs["embedding_keys"] is None
    assert kwargs["knn_key"] is None
    assert kwargs["gene_id_column"] is None
    assert kwargs["expression_source"] == "X"


def _prepared_export(tmp_path: Path) -> Path:
    export = tmp_path / "export"
    export.mkdir()
    (export / "dataset_identity.json").write_text(
        json.dumps({"version": 3, "id": "e", "name": "E"}),
        encoding="utf-8",
    )
    (export / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (export / "embeddings").mkdir()
    (export / "embeddings" / "0_2d.bin").write_bytes(b"\x00" * 8)
    return export


@pytest.mark.parametrize(
    ("flag", "value"),
    [
        ("--embedding-1d", "X_line"),
        ("--embedding-2d", "X_pca:0,1"),
        ("--embedding-3d", "X_umap_3d"),
        ("--knn", "connectivities"),
        ("--gene-id-column", "gene_symbols"),
        ("--expression", "raw"),
    ],
)
def test_every_mapping_flag_is_refused_for_a_prepared_export(
    tmp_path: Path,
    flag: str,
    value: str,
) -> None:
    args = create_parser().parse_args(["serve", str(_prepared_export(tmp_path)), "--quiet", flag, value])
    with (
        mock.patch("cellucid.server.serve") as serve,
        pytest.raises(ValueError, match=f"{flag} may only be used with direct AnnData input"),
    ):
        _run_serve(args)
    serve.assert_not_called()


def test_the_replaced_connectivity_switch_no_longer_parses() -> None:
    with pytest.raises(SystemExit):
        create_parser().parse_args(["serve", "data.h5ad", "--connectivity"])
