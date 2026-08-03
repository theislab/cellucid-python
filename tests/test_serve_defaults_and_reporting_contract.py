"""What ``cellucid serve`` decides on its own, and what it says while deciding.

Three rules meet here. An embedding is read from the key that names its
dimension, or from Scanpy's own ``X_umap`` at the dimension its column count
states. A neighbor graph is read only when it is asked for, because reading one
is the longest thing opening a large object can do. And a bind to every
interface is never rendered as a URL, because no client can open one.
"""

from __future__ import annotations

import ipaddress
import socket
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from cellucid._server_base import (
    bind_host_is_wildcard,
    format_http_origin,
    loopback_origin_for_bind,
    machine_origin,
)
from cellucid.anndata_adapter import AnnDataAdapter
from cellucid.anndata_server import AnnDataServer
from cellucid.connectivity_contract import validate_connectivity_edges

N_CELLS = 8


def _adata(*, obsm: dict[str, Any] | None = None, connectivity: bool = False) -> ad.AnnData:
    rng = np.random.default_rng(11)
    adata = ad.AnnData(
        X=rng.random((N_CELLS, 2), dtype=np.float32),
        obs=pd.DataFrame(
            {"leiden": pd.Categorical(["0", "1"] * (N_CELLS // 2))},
            index=[f"cell-{index}" for index in range(N_CELLS)],
        ),
        var=pd.DataFrame(index=["GeneA", "GeneB"]),
        obsm=obsm
        if obsm is not None
        else {
            "X_umap": rng.random((N_CELLS, 2)).astype(np.float32),
            "velocity_umap_2d": rng.random((N_CELLS, 2)).astype(np.float32),
        },
    )
    if connectivity:
        neighbors = np.zeros((N_CELLS, N_CELLS), dtype=np.float64)
        for index in range(N_CELLS - 1):
            neighbors[index, index + 1] = 0.5
            neighbors[index + 1, index] = 0.5
        adata.obsp["connectivities"] = sparse.csr_matrix(neighbors)
    return adata


def _adapter(adata: ad.AnnData, **kwargs: Any) -> AnnDataAdapter:
    return AnnDataAdapter(
        adata,
        dataset_name="Serve defaults",
        dataset_id="serve-defaults",
        centroid_min_points=1,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# A graph is read only when it is asked for
# ---------------------------------------------------------------------------


def test_a_graph_is_not_read_when_it_was_not_asked_for() -> None:
    adapter = _adapter(_adata(connectivity=True))
    try:
        assert adapter.serve_connectivity is False
        assert adapter.has_connectivity() is False
        assert adapter.get_connectivity_manifest() is None
        identity = adapter.get_dataset_identity()
        # The browser refuses a dataset whose flag and manifest disagree, and
        # refuses a zero it reads as a count, so this must be null and not 0.
        assert identity["stats"]["has_connectivity"] is False
        assert identity["stats"]["n_edges"] is None
    finally:
        adapter.close()


def test_the_same_graph_is_served_in_full_when_it_is_asked_for() -> None:
    adapter = _adapter(_adata(connectivity=True), serve_connectivity=True)
    try:
        assert adapter.has_connectivity() is True
        manifest = adapter.get_connectivity_manifest()
        assert manifest is not None
        assert manifest["n_edges"] == N_CELLS - 1
        assert adapter.get_dataset_identity()["stats"]["n_edges"] == N_CELLS - 1
    finally:
        adapter.close()


def test_asking_for_a_graph_that_does_not_exist_is_an_error() -> None:
    with pytest.raises(ValueError, match="Connectivity was asked for"):
        _adapter(_adata(connectivity=False), serve_connectivity=True)


def test_an_unasked_graph_is_never_touched_by_construction() -> None:
    """Opting out has to mean not reading obsp, not reading it and discarding it."""
    reads: list[str] = []

    class _CountingObsp(dict):
        def __getitem__(self, key: str) -> Any:
            reads.append(key)
            return super().__getitem__(key)

    adata = _adata(connectivity=True)
    graph = adata.obsp["connectivities"]
    object.__setattr__(adata, "_obsp", _CountingObsp({"connectivities": graph}))

    adapter = _adapter(adata)
    try:
        assert reads == []
    finally:
        adapter.close()


def test_the_dataset_report_separates_no_graph_from_an_unserved_one(capsys) -> None:
    server = AnnDataServer(
        _adata(connectivity=True),
        host="127.0.0.1",
        port=0,
        quiet=False,
        serve_web_ui=False,
        dataset_name="Report",
        dataset_id="report",
    )
    try:
        printed = capsys.readouterr().out
        assert "Connectivity: not served" in printed
        assert "--connectivity" in printed
    finally:
        server.stop()

    server = AnnDataServer(
        _adata(connectivity=False),
        host="127.0.0.1",
        port=0,
        quiet=False,
        serve_web_ui=False,
        dataset_name="Report",
        dataset_id="report",
    )
    try:
        assert "Connectivity: no\n" in capsys.readouterr().out
    finally:
        server.stop()


# ---------------------------------------------------------------------------
# Degree counting
# ---------------------------------------------------------------------------


def test_degrees_count_every_cell_including_one_no_edge_reaches() -> None:
    """``minlength`` is what keeps an isolated last cell in the degree array."""
    graph = np.zeros((6, 6), dtype=np.float64)
    graph[0, 1] = graph[1, 0] = 1.0
    edges = validate_connectivity_edges(sparse.csr_matrix(graph), n_cells=6)

    assert edges.n_edges == 1
    assert edges.max_neighbors == 1
    # A native int, because the manifest contract accepts nothing else.
    assert type(edges.max_neighbors) is int


def test_the_neighbor_maximum_is_the_true_maximum_degree() -> None:
    graph = np.zeros((5, 5), dtype=np.float64)
    for other in (1, 2, 3):
        graph[0, other] = graph[other, 0] = 1.0
    edges = validate_connectivity_edges(sparse.csr_matrix(graph), n_cells=5)

    assert edges.n_edges == 3
    assert edges.max_neighbors == 3


def test_a_large_graph_names_its_own_cost_before_paying_it(caplog, monkeypatch) -> None:
    """The size is reported before the work its size implies, not after it."""
    from cellucid import connectivity_contract

    graph = np.zeros((4, 4), dtype=np.float64)
    graph[0, 1] = graph[1, 0] = 1.0
    graph[2, 3] = graph[3, 2] = 1.0
    matrix = sparse.csr_matrix(graph)

    monkeypatch.setattr(connectivity_contract, "_LARGE_GRAPH_STORED_NEIGHBORS", 1)
    with caplog.at_level("WARNING", logger="cellucid.connectivity"):
        edges = validate_connectivity_edges(matrix, n_cells=4)

    assert edges.n_edges == 2
    assert "stored neighbors" in caplog.text
    assert "Serve without connectivity" in caplog.text


def test_an_ordinary_graph_reports_nothing(caplog) -> None:
    graph = np.zeros((4, 4), dtype=np.float64)
    graph[0, 1] = graph[1, 0] = 1.0
    with caplog.at_level("WARNING", logger="cellucid.connectivity"):
        validate_connectivity_edges(sparse.csr_matrix(graph), n_cells=4)
    assert caplog.text == ""


# ---------------------------------------------------------------------------
# A bind is not a destination
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("host", ["0.0.0.0", "::", "0:0:0:0:0:0:0:0", ""])
def test_every_wildcard_bind_is_recognized(host: str) -> None:
    assert bind_host_is_wildcard(host) is True


@pytest.mark.parametrize("host", ["127.0.0.1", "::1", "10.0.0.5", "localhost", "example.org"])
def test_no_reachable_host_is_mistaken_for_a_wildcard(host: str) -> None:
    assert bind_host_is_wildcard(host) is False


def test_an_ipv6_origin_is_bracketed_as_a_url_requires() -> None:
    assert format_http_origin("::1", 8765) == "http://[::1]:8765"
    assert format_http_origin("127.0.0.1", 8765) == "http://127.0.0.1:8765"
    assert format_http_origin("localhost", 8765) == "http://localhost:8765"


def test_a_wildcard_bind_reports_the_loopback_of_its_own_family() -> None:
    assert loopback_origin_for_bind("0.0.0.0", 8765) == "http://127.0.0.1:8765"
    assert loopback_origin_for_bind("", 8765) == "http://127.0.0.1:8765"
    assert loopback_origin_for_bind("::", 8765) == "http://[::1]:8765"


def test_a_machine_origin_is_a_real_name_or_nothing() -> None:
    origin = machine_origin(8765)
    if origin is None:
        return
    assert origin.startswith("http://")
    assert origin.endswith(":8765")
    host = origin[len("http://") : -len(":8765")].strip("[]")
    try:
        assert not ipaddress.ip_address(host).is_loopback
    except ValueError:
        assert host != "localhost"


def test_a_wildcard_bound_server_never_publishes_the_wildcard_as_a_url() -> None:
    server = AnnDataServer(
        _adata(),
        host="0.0.0.0",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Wildcard",
        dataset_id="wildcard",
    )
    server.start_background()
    try:
        assert "0.0.0.0" not in server.url
        assert "0.0.0.0" not in server.viewer_url
        assert server.url == f"http://127.0.0.1:{server.port}"
        assert server.viewer_url.endswith("/?anndata=true")
        for label, value in server.network_urls:
            assert "0.0.0.0" not in value
            assert label.strip().endswith(":")
        # The bind address itself is still reported as what it is.
        assert server.server_info["host"] == "0.0.0.0"
    finally:
        server.stop()


def test_a_loopback_bound_server_reports_no_network_urls() -> None:
    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Loopback",
        dataset_id="loopback",
    )
    server.start_background()
    try:
        assert server.network_urls == []
        assert server.url == f"http://127.0.0.1:{server.port}"
    finally:
        server.stop()


def test_the_banner_names_the_wildcard_in_words_and_never_in_a_url(capsys) -> None:
    from cellucid._server_base import print_server_banner

    print_server_banner(
        "http://127.0.0.1:8765",
        "http://127.0.0.1:8765/?anndata=true",
        network_urls=[
            ("Machine URL:  ", f"http://{socket.gethostname()}:8765"),
        ],
    )
    printed = capsys.readouterr().out
    assert "Bound to every network interface" in printed
    assert "0.0.0.0" not in printed
    assert "Local URL:    http://127.0.0.1:8765" in printed


def test_a_banner_without_network_urls_is_byte_identical_to_the_loopback_form(capsys) -> None:
    from cellucid._server_base import print_server_banner

    print_server_banner("http://127.0.0.1:8765", "http://127.0.0.1:8765/?anndata=true")
    printed = capsys.readouterr().out
    assert "Bound to every network interface" not in printed
    assert printed.count("Viewer URL:") == 1


# ---------------------------------------------------------------------------
# What the terminal says while a long build runs
# ---------------------------------------------------------------------------


def test_opening_an_object_reports_every_phase_of_the_build(capsys) -> None:
    """A build that can run for minutes must not run silently.

    Each line here marks a phase that is unbounded on a large object: the obs
    classification pass, the embedding resolution, the vector-field scan, and
    the manifest and centroid pass that used to run after a success line with
    nothing printed at all.
    """
    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=False,
        serve_web_ui=False,
        dataset_name="Progress",
        dataset_id="progress",
    )
    try:
        printed = capsys.readouterr().out
    finally:
        server.stop()

    for expected in [
        "[1/5] Detecting format",
        "[2/5] Loading AnnData",
        "Obs columns: classifying",
        "Embeddings: resolving obsm keys",
        "[3/5] Analyzing dataset",
        "[4/5] Building manifests",
        "Centroids: one per categorical field and dimension",
    ]:
        assert expected in printed, expected

    # The resolved key is named, so a reader can see which array was drawn.
    assert "Embeddings: 2D from obsm['X_umap']" in printed
    # Every step numbers itself out of the same total.
    assert "/4]" not in printed


def test_a_quiet_build_reports_nothing_at_all(capsys) -> None:
    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Quiet",
        dataset_id="quiet",
    )
    try:
        assert capsys.readouterr().out == ""
    finally:
        server.stop()


def test_the_adapter_is_silent_unless_a_caller_asks_for_progress(capsys) -> None:
    """Used as a library, the adapter prints nothing by default."""
    adapter = _adapter(_adata())
    try:
        assert adapter.quiet is True
        assert capsys.readouterr().out == ""
    finally:
        adapter.close()

    adapter = _adapter(_adata(), quiet=False)
    try:
        assert "Embeddings" in capsys.readouterr().out
    finally:
        adapter.close()


def test_asking_for_a_graph_reports_it_before_reading_it(capsys) -> None:
    adapter = _adapter(_adata(connectivity=True), serve_connectivity=True, quiet=False)
    try:
        printed = capsys.readouterr().out
        assert "Connectivity: reading obsp['connectivities']" in printed
        assert "edges" in printed
    finally:
        adapter.close()


# ---------------------------------------------------------------------------
# The CLI flag
# ---------------------------------------------------------------------------


def test_the_connectivity_flag_is_off_unless_it_is_given() -> None:
    from cellucid.cli import create_parser

    base = ["serve", "data.h5ad", "--dataset-name", "N", "--dataset-id", "n"]
    assert create_parser().parse_args(base).connectivity is None
    assert create_parser().parse_args([*base, "--connectivity"]).connectivity is True


def test_the_connectivity_flag_is_refused_for_a_prepared_export(tmp_path) -> None:
    import json

    from cellucid.cli import _run_serve, create_parser

    export = tmp_path / "export"
    export.mkdir()
    (export / "dataset_identity.json").write_text(
        json.dumps({"version": 2, "id": "e", "name": "E"}),
        encoding="utf-8",
    )
    (export / "obs_manifest.json").write_text("{}", encoding="utf-8")
    (export / "points_2d.bin").write_bytes(b"\x00" * 8)

    args = create_parser().parse_args(["serve", str(export), "--quiet", "--connectivity"])
    with pytest.raises(ValueError, match="--connectivity"):
        _run_serve(args)


# ---------------------------------------------------------------------------
# Vector fields are asked for on the same terms as connectivity
# ---------------------------------------------------------------------------


def test_vector_fields_are_not_read_when_they_were_not_asked_for() -> None:
    adapter = _adapter(_adata())
    try:
        assert adapter.serve_vector_fields is False
        assert "vector_fields" not in adapter.get_dataset_identity()
    finally:
        adapter.close()


def test_the_same_vector_fields_are_served_when_they_are_asked_for() -> None:
    adapter = _adapter(_adata(), serve_vector_fields=True)
    try:
        identity = adapter.get_dataset_identity()
        assert list(identity["vector_fields"]["fields"]) == ["velocity_umap"]
    finally:
        adapter.close()


def test_asking_for_vector_fields_an_object_declares_none_of_is_not_an_error() -> None:
    """A key outside the grammar is not a malformed field: it is not a field.

    This is where vector fields part from connectivity, which names one exact
    matrix and therefore cannot be absent once it has been asked for.
    """
    adata = _adata(obsm={"X_umap": np.zeros((N_CELLS, 2), dtype=np.float32) + np.arange(
        N_CELLS, dtype=np.float32
    )[:, None]})
    adapter = _adapter(adata, serve_vector_fields=True)
    try:
        assert "vector_fields" not in adapter.get_dataset_identity()
    finally:
        adapter.close()


def test_a_default_cannot_name_fields_nobody_asked_to_serve() -> None:
    with pytest.raises(ValueError, match="vector_field_default"):
        _adapter(_adata(), vector_field_default="velocity_umap")


def test_the_report_separates_no_vector_fields_from_unserved_ones(capsys) -> None:
    server = AnnDataServer(
        _adata(),
        host="127.0.0.1",
        port=0,
        quiet=False,
        serve_web_ui=False,
        dataset_name="Report",
        dataset_id="report",
    )
    try:
        printed = capsys.readouterr().out
        assert "Vector fields: not served" in printed
        assert "--vector-fields" in printed
    finally:
        server.stop()


def test_the_vector_fields_flag_is_off_unless_it_is_given() -> None:
    from cellucid.cli import create_parser

    base = ["serve", "data.h5ad", "--dataset-name", "N", "--dataset-id", "n"]
    assert create_parser().parse_args(base).vector_fields is None
    assert create_parser().parse_args([*base, "--vector-fields"]).vector_fields is True
