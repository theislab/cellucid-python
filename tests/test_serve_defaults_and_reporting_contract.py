"""What ``cellucid serve`` decides on its own, and what it says while deciding.

Three rules meet here. An embedding is read from the key that names its
dimension, or from Scanpy's own ``X_umap`` at the dimension its column count
states, unless a slot is named explicitly. A neighbor graph is read only when
``knn_key`` names one, because reading one is the longest thing opening a large
object can do. And a bind to every
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
from cellucid.anndata_server import AnnDataServer

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
    """Every phase of a start is named, and an object says where it is written."""
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
        "[1/4] Detecting format",
        "Format: written once to a temporary .h5ad",
        "[2/4] Reading metadata",
        "obsm: X_umap",
        "[3/4] Checking launch choices",
    ]:
        assert expected in printed, expected
    # Every step numbers itself out of the same total.
    assert "/5]" not in printed


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


# ---------------------------------------------------------------------------
# The CLI flag
# ---------------------------------------------------------------------------


def test_the_knn_flag_is_off_unless_it_names_a_matrix() -> None:
    from cellucid.cli import create_parser

    base = ["serve", "data.h5ad", "--dataset-name", "N", "--dataset-id", "n"]
    assert create_parser().parse_args(base).knn is None
    assert create_parser().parse_args([*base, "--knn", "connectivities"]).knn == "connectivities"
    # The boolean flag it replaced is gone, not aliased.
    with pytest.raises(SystemExit):
        create_parser().parse_args([*base, "--connectivity"])


# ---------------------------------------------------------------------------
# Vector fields are asked for on the same terms as connectivity
# ---------------------------------------------------------------------------


def test_vector_fields_need_no_flag_to_be_drawn() -> None:
    """The viewer reads a vector field only when the overlay shows it."""
    from cellucid.cli import create_parser

    base = ["serve", "data.h5ad", "--dataset-name", "N", "--dataset-id", "n"]
    assert not hasattr(create_parser().parse_args(base), "vector_fields")
    with pytest.raises(SystemExit):
        create_parser().parse_args([*base, "--vector-fields"])
