"""The AnnData server's public signatures and its single-use lifecycle.

A server owns two resources: its socket, and -- for an in-memory AnnData
object -- the store in priority order it serves, with the files that store
lays a cell-major expression matrix out in when the viewer first reads it.
Every way a server ends (a stop, a failed start, a failed serve loop, a
constructor that fails, or Python dropping the object) releases both, and a
stopped server cannot be restarted.
"""

from __future__ import annotations

import gc
import inspect
import socket
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from _mapped_files import mapped_files
from anndata import AnnData
from scipy import sparse

from cellucid.anndata_server import AnnDataServer, serve_anndata
from cellucid.jupyter import AnnDataViewer, Widget, show_anndata


def _adata() -> AnnData:
    return AnnData(
        X=sparse.csr_matrix(np.array([[1.0], [2.0]], dtype=np.float32)),
        obs=pd.DataFrame(index=pd.Index(["cell-1", "cell-2"], dtype=object)),
        var=pd.DataFrame(index=pd.Index(["gene-1"], dtype=object)),
        obsm={"X_umap_2d": np.array([[0.0, 0.0], [1.0, 1.0]], dtype=np.float32)},
    )


def _server(data: object = None, **options: object) -> AnnDataServer:
    return AnnDataServer(
        _adata() if data is None else data,
        host="127.0.0.1",
        port=options.pop("port", 0),
        quiet=True,
        serve_web_ui=options.pop("serve_web_ui", False),
        dataset_name="Lifecycle fixture",
        dataset_id=options.pop("dataset_id", "lifecycle-fixture"),
        **options,
    )


def _assert_socket_released(port: int) -> None:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", port))


def _held(server: AnnDataServer) -> Path:
    """The folder the served store laid the object's cell-major ``X`` out in, made by reading it."""
    store = server._served.priority
    assert store is not None
    assert store.get("X/indptr/0") is not None
    folder = store._folder.folder
    assert folder is not None and folder.is_dir()
    return folder


def _assert_released(server: AnnDataServer, temporary: Path) -> None:
    # Every mapping let go first: Windows removes no mapped file.
    assert mapped_files(temporary) == []
    assert not temporary.exists()
    assert server._server is None
    assert server.is_running() is False
    with pytest.raises(RuntimeError, match="single-use"):
        server.start_background()


LAUNCH_CHOICES = [
    "latent_key",
    "gene_id_column",
    "normalize_embeddings",
    "centroid_min_points",
    "dataset_name",
    "dataset_id",
    "obs_keys",
    "embedding_keys",
    "knn_key",
    "expression_source",
    "vector_field_default",
    "modality",
]


def test_public_anndata_entry_points_have_complete_closed_signatures() -> None:
    server_parameters = [
        "data",
        "port",
        "host",
        "open_browser",
        "quiet",
        *LAUNCH_CHOICES,
        "serve_web_ui",
        "web_source_url",
        "web_cache_dir",
        "allowed_hosts",
    ]
    assert list(inspect.signature(AnnDataServer).parameters) == server_parameters
    assert list(inspect.signature(serve_anndata).parameters) == server_parameters
    notebook_tail = ["client_server_url", "web_source_url", "web_cache_dir", "allowed_hosts"]
    assert list(inspect.signature(AnnDataViewer).parameters) == [
        "data",
        "port",
        "auto_open",
        *LAUNCH_CHOICES,
        *notebook_tail,
    ]
    # The notebook widget runs no server: the browser loads the app from web_source_url.
    widget_parameters = ["data", "height", *LAUNCH_CHOICES, "web_source_url"]
    assert list(inspect.signature(show_anndata).parameters) == widget_parameters
    assert list(inspect.signature(Widget).parameters) == widget_parameters


def test_show_anndata_rejects_identity_before_viewer_setup() -> None:
    with (
        mock.patch(
            "cellucid.jupyter.AnnDataViewer._init_viewer_state",
            side_effect=AssertionError("Notebook viewer setup was entered"),
        ) as viewer_setup,
        pytest.raises(ValueError, match="dataset_id"),
    ):
        show_anndata(object(), dataset_name="Declared name", dataset_id="data/set")

    viewer_setup.assert_not_called()


def test_a_file_on_disk_is_read_where_it_is_and_its_decoded_arrays_go_with_the_server(tmp_path: Path) -> None:
    source = tmp_path / "fixture.h5ad"
    _adata().write_h5ad(source)
    server = _server(source)
    store = server._served.priority
    assert store is not None
    # A cell-major X on disk is served cell-major, read from the file: an array
    # stored in chunks (anndata's sparse arrays) is decoded once into the
    # scratch folder as its rows are read, an array in one piece read where it
    # lies; nothing else is written, and the folder goes with the server.
    assert store.get("X/indptr/0") is not None and store.get("X/data/0") is not None
    folder = store._folder.folder
    assert folder is not None
    assert {path.name.split("-")[0] for path in folder.iterdir()} <= {"owner", "source"}
    server.stop()
    assert mapped_files(folder) == [] and mapped_files(source) == []
    assert not folder.exists()
    assert source.exists()


def test_an_object_is_served_from_a_store_released_on_stop() -> None:
    server = _server()
    temporary = _held(server)
    server.start_background()
    port = server.port
    server.stop()

    _assert_released(server, temporary)
    _assert_socket_released(port)
    with pytest.raises(RuntimeError, match="not running"):
        _ = server.url
    with pytest.raises(RuntimeError, match="not running"):
        _ = server.viewer_url


def test_stop_before_start_consumes_the_server_and_its_file() -> None:
    server = _server()
    temporary = _held(server)
    server.stop()
    _assert_released(server, temporary)


def test_a_server_python_drops_without_a_stop_removes_its_file() -> None:
    server = _server()
    temporary = _held(server)
    del server
    gc.collect()
    assert mapped_files(temporary) == []
    assert not temporary.exists()


def test_a_constructor_that_fails_releases_the_store_it_made() -> None:
    from cellucid._served_store import PriorityStore

    made: list[PriorityStore] = []
    real_of_anndata = PriorityStore.of_anndata

    def tracking(adata: AnnData) -> PriorityStore:
        store = real_of_anndata(adata)
        made.append(store)
        assert store.get("X/indptr/0") is not None
        return store

    with (
        mock.patch.object(PriorityStore, "of_anndata", tracking),
        pytest.raises(ValueError, match="knn_key names obsp"),
    ):
        _server(knn_key="connectivities")
    assert len(made) == 1
    assert made[0]._folder.folder is not None and mapped_files(made[0]._folder.folder) == []
    assert not made[0]._folder.folder.exists()
    assert not made[0]._finalizer.alive


def test_bind_failure_releases_the_file_and_consumes_the_server() -> None:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as occupied:
        occupied.bind(("127.0.0.1", 0))
        occupied.listen(1)
        server = _server(port=occupied.getsockname()[1])
        temporary = _held(server)
        with pytest.raises(OSError):
            server.start_background()
    _assert_released(server, temporary)


def test_web_cache_failure_preserves_its_error_and_releases_the_file(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = _server(serve_web_ui=True)
    temporary = _held(server)
    original = RuntimeError("exact cache failure")

    def fail_cache(**kwargs: object) -> None:
        assert kwargs["force"] is True
        raise original

    monkeypatch.setattr("cellucid.web_cache.ensure_web_ui_cached", fail_cache)
    with pytest.raises(RuntimeError, match="exact cache failure") as raised:
        server.start_background()

    assert raised.value is original
    _assert_released(server, temporary)


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
        assert ("Viewer UI generation", "establishing exact configured source") in detail_calls
        assert "Viewer UI generation established" in success_calls
        assert len(cache_calls) == 1
        assert cache_calls[0]["force"] is True
        assert cache_calls[0]["source_url"] == "https://viewer.example"
    finally:
        server.stop()


def test_thread_start_failure_rolls_back_the_socket_and_the_file(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = _server()
    temporary = _held(server)
    original = RuntimeError("exact thread start failure")

    def fail_start(_thread: object) -> None:
        raise original

    monkeypatch.setattr("cellucid.anndata_server.threading.Thread.start", fail_start)
    with pytest.raises(RuntimeError, match="exact thread start failure") as raised:
        server.start_background()

    assert raised.value is original
    _assert_released(server, temporary)
    _assert_socket_released(server.port)


@pytest.mark.parametrize("browser", ["false", "raises"])
def test_a_browser_that_does_not_open_is_a_transactional_start_failure(
    browser: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = _server(open_browser=True)
    temporary = _held(server)
    original = RuntimeError("exact browser failure")

    def open_browser(_url: str) -> bool:
        if browser == "raises":
            raise original
        return False

    monkeypatch.setattr("cellucid.anndata_server.webbrowser.open", open_browser)
    with pytest.raises(RuntimeError, match="exact browser failure|Could not open the browser"):
        server.start_background()

    _assert_released(server, temporary)
    _assert_socket_released(server.port)


def test_background_serve_failure_is_retained_for_wait_and_releases_resources(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    server = _server()
    temporary = _held(server)
    original = RuntimeError("exact serve failure")

    def fail_serve(_server: object) -> None:
        raise original

    monkeypatch.setattr(
        "cellucid.anndata_server.CellucidThreadingHTTPServer.serve_forever",
        fail_serve,
    )
    server.start_background()
    with pytest.raises(RuntimeError, match="exact serve failure") as raised:
        server.wait()

    assert raised.value is original
    assert mapped_files(temporary) == []
    assert not temporary.exists()
    assert server._server is None
    assert server.is_running() is False
    _assert_socket_released(server.port)


def test_stopping_closes_every_kept_alive_connection() -> None:
    """A browser's kept connection must not outlive the server that answered it.

    Closing only the listening socket left each kept connection answered by
    its handler thread: a server started next on the same port had the page
    read the stopped server's dataset documents over those connections.
    """
    import http.client
    import json

    first = _server(dataset_id="first-dataset")
    first.start_background()
    port = first.port
    kept = http.client.HTTPConnection("127.0.0.1", port, timeout=10)
    kept.request("GET", "/_cellucid/datasets")
    assert json.loads(kept.getresponse().read())["datasets"][0]["id"] == "first-dataset"
    first.stop()

    second = _server(port=port, dataset_id="second-dataset")
    second.start_background()
    try:
        try:
            kept.request("GET", "/_cellucid/datasets")
            answered = json.loads(kept.getresponse().read())["datasets"][0]["id"]
        except (ConnectionError, http.client.HTTPException):
            answered = None
        assert answered != "first-dataset"
    finally:
        kept.close()
        second.stop()


def test_a_notebook_viewer_takes_the_next_free_port_when_its_own_is_taken() -> None:
    """As cellucid-r's cellucid_view(): 8765 (Google sign-in's port), else the next free one."""
    holder = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    holder.bind(("127.0.0.1", 0))
    holder.listen(1)
    taken = holder.getsockname()[1]
    try:
        refused = _server(port=taken)
        with pytest.raises(OSError):
            refused.start_background()
        walking = _server(port=taken)
        walking.next_free_port = True
        walking.start_background()
        try:
            assert walking.port > taken
        finally:
            walking.stop()
    finally:
        holder.close()
    viewer_port = inspect.signature(AnnDataViewer.__init__).parameters["port"].default
    assert viewer_port is None  # None is 8765 with the walk, a given port is used or refused
