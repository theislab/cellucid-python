from __future__ import annotations

import inspect
import socket
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.anndata_server import AnnDataServer, serve_anndata
from cellucid.jupyter import AnnDataViewer, show_anndata


def _adata() -> AnnData:
    return AnnData(
        X=np.array([[1.0], [2.0]], dtype=np.float32),
        obs=pd.DataFrame(
            index=pd.Index(["cell-1", "cell-2"], dtype=object)
        ),
        var=pd.DataFrame(index=pd.Index(["gene-1"], dtype=object)),
        obsm={
            "X_umap_2d": np.array(
                [[0.0, 0.0], [1.0, 1.0]],
                dtype=np.float32,
            )
        },
    )


def _write_h5ad(path: Path) -> None:
    _adata().write_h5ad(path)


def _backed_server(path: Path) -> AnnDataServer:
    return AnnDataServer(
        path,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Lifecycle fixture",
        dataset_id="lifecycle-fixture",
    )


def _assert_socket_released(port: int) -> None:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", port))


def test_public_anndata_entry_points_have_complete_closed_signatures() -> None:
    assert list(inspect.signature(AnnDataAdapter.from_file).parameters) == [
        "path",
        "latent_key",
        "gene_id_column",
        "normalize_embeddings",
        "centroid_outlier_quantile",
        "centroid_min_points",
        "dataset_name",
        "dataset_id",
        "obs_keys",
        "vector_field_default",
        "serve_connectivity",
        "serve_vector_fields",
        "quiet",
    ]
    server_parameters = [
        "data",
        "port",
        "host",
        "open_browser",
        "quiet",
        "latent_key",
        "gene_id_column",
        "normalize_embeddings",
        "centroid_outlier_quantile",
        "centroid_min_points",
        "dataset_name",
        "dataset_id",
        "obs_keys",
        "vector_field_default",
        "serve_connectivity",
        "serve_vector_fields",
        "serve_web_ui",
        "web_source_url",
        "web_cache_dir",
        "allowed_hosts",
    ]
    assert list(inspect.signature(AnnDataServer).parameters) == server_parameters
    assert list(inspect.signature(serve_anndata).parameters) == server_parameters
    assert list(inspect.signature(AnnDataViewer).parameters) == [
        "data",
        "port",
        "height",
        "auto_open",
        "latent_key",
        "gene_id_column",
        "normalize_embeddings",
        "centroid_outlier_quantile",
        "centroid_min_points",
        "dataset_name",
        "dataset_id",
        "obs_keys",
        "vector_field_default",
        "serve_connectivity",
        "serve_vector_fields",
        "client_server_url",
        "web_source_url",
        "web_cache_dir",
        "allowed_hosts",
    ]
    assert list(inspect.signature(show_anndata).parameters) == [
        "data",
        "height",
        "latent_key",
        "gene_id_column",
        "normalize_embeddings",
        "centroid_outlier_quantile",
        "centroid_min_points",
        "dataset_name",
        "dataset_id",
        "obs_keys",
        "vector_field_default",
        "serve_connectivity",
        "serve_vector_fields",
        "client_server_url",
        "web_source_url",
        "web_cache_dir",
        "allowed_hosts",
    ]


def test_show_anndata_rejects_identity_before_base_viewer_setup() -> None:
    with (
        mock.patch(
            "cellucid.jupyter.BaseViewer.__init__",
            side_effect=AssertionError("Notebook viewer setup was entered"),
        ) as base_init,
        pytest.raises(ValueError, match="dataset_id"),
    ):
        show_anndata(
            object(),
            dataset_name="Declared name",
            dataset_id="data/set",
        )

    base_init.assert_not_called()


@pytest.mark.parametrize("backed_value", [None, True, False, "r", "r+"])
def test_removed_backed_modes_are_not_accepted(
    tmp_path: Path,
    backed_value: object,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)

    with pytest.raises(TypeError, match="backed"):
        AnnDataAdapter.from_file(
            source,
            dataset_name="Exact read owner",
            dataset_id="exact-read-owner",
            backed=backed_value,
        )


def test_h5ad_file_loading_has_one_read_only_backed_owner(tmp_path: Path) -> None:
    source = tmp_path / "fixture.h5ad"
    source.write_bytes(b"reader dispatch is mocked")
    adata = _adata()

    with mock.patch("anndata.read_h5ad", return_value=adata) as read_h5ad:
        adapter = AnnDataAdapter.from_file(
            source,
            dataset_name="Exact read owner",
            dataset_id="exact-read-owner",
        )
    try:
        read_h5ad.assert_called_once_with(source.resolve(), backed="r")
    finally:
        adapter.close()


def test_write_enabled_backed_anndata_is_rejected(tmp_path: Path) -> None:
    import anndata

    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)
    write_enabled = anndata.read_h5ad(source, backed="r+")
    try:
        with pytest.raises(ValueError, match="opened read-only"):
            AnnDataAdapter(
                write_enabled,
                dataset_name="Write-enabled input",
                dataset_id="write-enabled-input",
            )
        assert write_enabled.file.is_open is True
    finally:
        write_enabled.file.close()


def test_stopped_server_is_single_use_and_capabilities_reject_after_close() -> None:
    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Single use",
        dataset_id="single-use",
    )
    server.start_background()
    port = server.port
    server.stop()

    assert server.is_running() is False
    assert server._server is None
    _assert_socket_released(port)
    with pytest.raises(RuntimeError, match="closed"):
        _ = server.adapter.is_backed
    with pytest.raises(RuntimeError, match="closed"):
        server.adapter.has_connectivity()
    with pytest.raises(RuntimeError, match="not running"):
        _ = server.url
    with pytest.raises(RuntimeError, match="not running"):
        _ = server.viewer_url
    with pytest.raises(RuntimeError, match="single-use"):
        server.start_background()


def test_stop_before_start_consumes_and_closes_the_server() -> None:
    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Stopped before start",
        dataset_id="stopped-before-start",
    )
    server.stop()

    with pytest.raises(RuntimeError, match="single-use"):
        server.start_background()
    with pytest.raises(RuntimeError, match="closed"):
        _ = server.adapter.is_backed


def test_bind_failure_closes_owned_backed_file_and_consumes_server(
    tmp_path: Path,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)

    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as occupied:
        occupied.bind(("127.0.0.1", 0))
        occupied.listen(1)
        port = occupied.getsockname()[1]
        server = AnnDataServer(
            source,
            host="127.0.0.1",
            port=port,
            quiet=True,
            serve_web_ui=False,
            dataset_name="Bind failure",
            dataset_id="bind-failure",
        )
        owned_adata = server.adapter.adata
        assert owned_adata.file.is_open is True

        with pytest.raises(OSError) as raised:
            server.start_background()

    assert raised.value is not None
    assert owned_adata.file.is_open is False
    assert server._server is None
    assert server.is_running() is False
    with pytest.raises(RuntimeError, match="single-use"):
        server.start_background()


def test_web_cache_failure_preserves_error_and_closes_backed_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)
    server = AnnDataServer(
        source,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=True,
        dataset_name="Cache failure",
        dataset_id="cache-failure",
    )
    owned_adata = server.adapter.adata
    original = RuntimeError("exact cache failure")

    def fail_cache(**kwargs: object) -> None:
        assert kwargs["force"] is True
        raise original

    monkeypatch.setattr("cellucid.web_cache.ensure_web_ui_cached", fail_cache)
    with pytest.raises(RuntimeError, match="exact cache failure") as raised:
        server.start_background()

    assert raised.value is original
    assert owned_adata.file.is_open is False
    assert server._server is None
    assert server.is_running() is False


def test_anndata_server_reports_exact_viewer_generation_establishment(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detail_calls: list[tuple[str, str]] = []
    success_calls: list[str] = []
    cache_calls: list[dict[str, object]] = []

    monkeypatch.setattr(
        "cellucid.anndata_server.print_detail",
        lambda label, value: detail_calls.append((label, value)),
    )
    monkeypatch.setattr(
        "cellucid.anndata_server.print_success",
        lambda message: success_calls.append(message),
    )
    monkeypatch.setattr(
        "cellucid.web_cache.ensure_web_ui_cached",
        lambda **kwargs: cache_calls.append(kwargs),
    )
    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=False,
        serve_web_ui=True,
        web_source_url="https://viewer.example",
        web_cache_dir=tmp_path / "web-cache",
        dataset_name="Generation report",
        dataset_id="generation-report",
    )

    server.start_background()
    try:
        assert (
            "Viewer UI generation",
            "establishing exact configured source",
        ) in detail_calls
        assert "Viewer UI generation established" in success_calls
        assert len(cache_calls) == 1
        assert cache_calls[0]["force"] is True
        assert cache_calls[0]["source_url"] == "https://viewer.example"
    finally:
        server.stop()


def test_constructor_contract_failure_preserves_error_and_closes_backed_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)
    import anndata

    real_read_h5ad = anndata.read_h5ad
    owned: list[AnnData] = []
    original = RuntimeError("exact contract construction failure")

    def tracking_read_h5ad(*args: object, **kwargs: object) -> AnnData:
        adata = real_read_h5ad(*args, **kwargs)
        owned.append(adata)
        return adata

    def fail_contract(
        _server: AnnDataServer,
    ) -> tuple[dict, dict[str, bytes], dict[str, int]]:
        raise original

    monkeypatch.setattr(anndata, "read_h5ad", tracking_read_h5ad)
    monkeypatch.setattr(AnnDataServer, "_build_http_contract", fail_contract)
    with pytest.raises(RuntimeError, match="exact contract construction failure") as raised:
        AnnDataServer(
            source,
            host="127.0.0.1",
            port=0,
            quiet=True,
            serve_web_ui=False,
            dataset_name="Constructor failure",
            dataset_id="constructor-failure",
        )

    assert raised.value is original
    assert len(owned) == 1
    assert owned[0].file.is_open is False


def test_thread_start_failure_rolls_back_bound_socket_and_backed_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)
    server = _backed_server(source)
    owned_adata = server.adapter.adata
    original = RuntimeError("exact thread start failure")

    def fail_start(_thread: object) -> None:
        raise original

    monkeypatch.setattr("cellucid.anndata_server.threading.Thread.start", fail_start)
    with pytest.raises(RuntimeError, match="exact thread start failure") as raised:
        server.start_background()

    assert raised.value is original
    assert owned_adata.file.is_open is False
    assert server._server is None
    assert server.is_running() is False
    _assert_socket_released(server.port)


def test_false_browser_result_is_a_transactional_start_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)
    server = AnnDataServer(
        source,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        open_browser=True,
        dataset_name="Browser result",
        dataset_id="browser-result",
    )
    owned_adata = server.adapter.adata
    monkeypatch.setattr(
        "cellucid.anndata_server.webbrowser.open",
        lambda _url: False,
    )

    with pytest.raises(RuntimeError, match="Could not open the browser"):
        server.start_background()

    assert owned_adata.file.is_open is False
    assert server._server is None
    assert server.is_running() is False
    _assert_socket_released(server.port)


def test_browser_failure_rolls_back_live_thread_socket_and_backed_file(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)
    server = AnnDataServer(
        source,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        open_browser=True,
        dataset_name="Browser failure",
        dataset_id="browser-failure",
    )
    owned_adata = server.adapter.adata
    original = RuntimeError("exact browser failure")

    def fail_browser(_url: str) -> bool:
        raise original

    monkeypatch.setattr("cellucid.anndata_server.webbrowser.open", fail_browser)
    with pytest.raises(RuntimeError, match="exact browser failure") as raised:
        server.start_background()

    assert raised.value is original
    assert owned_adata.file.is_open is False
    assert server._server is None
    assert server.is_running() is False
    _assert_socket_released(server.port)


def test_background_serve_failure_is_retained_for_wait_and_closes_resources(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "fixture.h5ad"
    _write_h5ad(source)
    server = _backed_server(source)
    owned_adata = server.adapter.adata
    original = RuntimeError("exact serve failure")

    def fail_serve(_server: object) -> None:
        raise original

    monkeypatch.setattr(
        "cellucid.anndata_server.HTTPServer.serve_forever",
        fail_serve,
    )
    server.start_background()
    with pytest.raises(RuntimeError, match="exact serve failure") as raised:
        server.wait()

    assert raised.value is original
    assert owned_adata.file.is_open is False
    assert server._server is None
    assert server.is_running() is False
    _assert_socket_released(server.port)
