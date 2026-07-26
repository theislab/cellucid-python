from __future__ import annotations

import http.client
import inspect
from collections.abc import Callable

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.anndata_server import AnnDataRequestHandler, AnnDataServer, serve_anndata
from cellucid.jupyter import AnnDataViewer, show_anndata


def _adata(*, index_column: list[str] | None = None) -> ad.AnnData:
    var = pd.DataFrame(index=["Index-A", "Index-B"])
    if index_column is not None:
        var["index"] = index_column
    adata = ad.AnnData(
        X=np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell-1", "cell-2"]),
        var=var,
    )
    adata.obsm["X_umap_2d"] = np.array(
        [[0.0, 0.0], [1.0, 1.0]],
        dtype=np.float32,
    )
    return adata


def _request(host: str, port: int, path: str) -> tuple[int, bytes]:
    connection = http.client.HTTPConnection(host, port, timeout=5)
    try:
        connection.request("GET", path)
        response = connection.getresponse()
        return response.status, response.read()
    finally:
        connection.close()


def test_direct_anndata_http_contract_has_no_dataset_prefixed_route_alias() -> None:
    assert "dataset_id" not in inspect.signature(AnnDataRequestHandler).parameters

    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Exact route",
        dataset_id="exact-route",
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        assert _request(host, port, "/dataset_identity.json")[0] == 200
        assert _request(host, port, "/exact-route/dataset_identity.json")[0] == 404
    finally:
        server.stop()


@pytest.mark.parametrize(
    ("entry_point", "parameter"),
    [
        (AnnDataAdapter, "gene_id_column"),
        (AnnDataAdapter.from_file, "gene_id_column"),
        (AnnDataServer, "gene_id_column"),
        (serve_anndata, "gene_id_column"),
        (AnnDataViewer, "gene_id_column"),
        (show_anndata, "gene_id_column"),
    ],
)
def test_direct_anndata_gene_index_selector_defaults_to_none(
    entry_point: Callable[..., object],
    parameter: str,
) -> None:
    assert inspect.signature(entry_point).parameters[parameter].default is None


def test_none_selects_var_index_and_is_declared_in_the_manifest() -> None:
    adapter = AnnDataAdapter(
        _adata(index_column=["Column-A", "Column-B"]),
        gene_id_column=None,
        dataset_name="Index contract",
        dataset_id="index-contract",
    )

    assert adapter.get_gene_ids() == ["Index-A", "Index-B"]
    assert adapter.get_var_manifest()["var_gene_id_column"] is None


def test_index_string_is_one_literal_var_column_name() -> None:
    adapter = AnnDataAdapter(
        _adata(index_column=["Column-A", "Column-B"]),
        gene_id_column="index",
        dataset_name="Literal column",
        dataset_id="literal-column",
    )

    assert adapter.get_gene_ids() == ["Column-A", "Column-B"]
    assert adapter.get_var_manifest()["var_gene_id_column"] == "index"


def test_index_string_is_rejected_when_that_exact_var_column_is_absent() -> None:
    with pytest.raises(
        KeyError,
        match=r"gene_id_column 'index' not found in var",
    ):
        AnnDataAdapter(
            _adata(),
            gene_id_column="index",
            dataset_name="Missing literal column",
            dataset_id="missing-literal-column",
        )


@pytest.mark.parametrize("gene_id_column", ["", "   "])
def test_gene_id_column_rejects_blank_strings(gene_id_column: str) -> None:
    with pytest.raises(ValueError, match="gene_id_column"):
        AnnDataAdapter(
            _adata(),
            gene_id_column=gene_id_column,
            dataset_name="Blank selector",
            dataset_id="blank-selector",
        )
