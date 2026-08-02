"""Both viewer servers publish the wire capabilities this installation accepts.

``GET /_cellucid/protocol`` exists so a web build can discover whether the
notebook driving it accepts an event before it emits one. The route therefore
has to answer identically from either server, has to answer without credentials
exactly as the other read-only ``/_cellucid/*`` routes do, and has to publish
the declaration :mod:`cellucid.jupyter._wire` validates against rather than a
copy kept in the server -- a published capability the validator would then
reject is worse than no endpoint at all.
"""

from __future__ import annotations

import http.client
import json
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import cellucid.jupyter._wire as wire
from cellucid.anndata_server import AnnDataServer
from cellucid.prepare_data import prepare
from cellucid.server import CellucidServer

PROTOCOL_PATH = "/_cellucid/protocol"


@contextmanager
def _prepared_server(tmp_path: Path) -> Iterator[tuple[str, int]]:
    output = tmp_path / "generation"
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.25], [0.25, 1.0]],
        dtype=np.float32,
    )
    prepare(
        latent_space=embedding.copy(),
        obs=pd.DataFrame({"score": [0.25, 0.5, 0.75]}),
        X_umap_2d=embedding,
        out_dir=output,
        dataset_name="Protocol endpoint contract",
        dataset_id="protocol-endpoint-contract",
        obs_categorical_dtype="uint16",
        centroid_min_points=1,
    )
    server = CellucidServer(
        output,
        host="127.0.0.1",
        port=0,
        open_browser=False,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        yield host, port
    finally:
        server.stop()


@contextmanager
def _anndata_server() -> Iterator[tuple[str, int]]:
    adata = ad.AnnData(
        X=np.array([[1.0], [2.0]], dtype=np.float32),
        obs=pd.DataFrame(index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["gene-1"]),
        obsm={
            "X_umap_2d": np.array(
                [[0.0, 0.0], [1.0, 1.0]],
                dtype=np.float32,
            )
        },
    )
    server = AnnDataServer(
        adata,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        dataset_name="Protocol endpoint contract",
        dataset_id="protocol-endpoint-contract",
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        yield host, port
    finally:
        server.stop()


@contextmanager
def _both_servers(tmp_path: Path) -> Iterator[dict[str, tuple[str, int]]]:
    with _prepared_server(tmp_path) as prepared, _anndata_server() as anndata:
        yield {"exported": prepared, "anndata": anndata}


def _request(
    address: tuple[str, int],
    path: str,
    *,
    method: str = "GET",
    headers: dict[str, str] | None = None,
) -> tuple[int, bytes, str | None, str | None]:
    """Send one exact request line without letting a client rewrite the path."""
    host, port = address
    connection = http.client.HTTPConnection(host, port, timeout=5)
    try:
        connection.request(method, path, headers=headers or {})
        response = connection.getresponse()
        body = response.read()
        return (
            response.status,
            body,
            response.getheader("Content-Type"),
            response.getheader("Content-Length"),
        )
    finally:
        connection.close()


@pytest.mark.parametrize("server_kind", ["exported", "anndata"])
def test_every_viewer_server_publishes_the_declaration_the_validator_reads(
    tmp_path: Path,
    server_kind: str,
) -> None:
    expected = wire._protocol_capability_document()
    with _both_servers(tmp_path) as servers:
        status, body, content_type, _length = _request(
            servers[server_kind],
            PROTOCOL_PATH,
        )

    assert status == 200
    assert content_type == "application/json"
    # Compared against the module the validator itself reads, never a literal:
    # a copy kept in the server is exactly the drift this route must not have.
    assert json.loads(body) == expected


def test_both_viewer_servers_publish_one_identical_capability_document(
    tmp_path: Path,
) -> None:
    with _both_servers(tmp_path) as servers:
        _status, exported_body, _type, _length = _request(
            servers["exported"],
            PROTOCOL_PATH,
        )
        _status, anndata_body, _type, _length = _request(
            servers["anndata"],
            PROTOCOL_PATH,
        )

    # A viewer asks this route before it knows which server answered, so the
    # two must not be distinguishable here.
    assert json.loads(exported_body) == json.loads(anndata_body)


@pytest.mark.parametrize("server_kind", ["exported", "anndata"])
def test_the_capability_document_carries_only_growable_capability_lists(
    tmp_path: Path,
    server_kind: str,
) -> None:
    with _both_servers(tmp_path) as servers:
        _status, body, _type, _length = _request(servers[server_kind], PROTOCOL_PATH)
    document: Any = json.loads(body)

    # The key set is fixed on purpose. A viewer validates payloads against exact
    # key sets, so a later capability has to arrive inside these lists rather
    # than as a new key, which would break the builds this route exists to
    # protect. There is deliberately no version to compare against either.
    assert set(document) == {"commands", "events"}
    for key, names in document.items():
        assert isinstance(names, list), key
        assert all(type(name) is str for name in names), key
        assert names == sorted(names), key
        assert len(set(names)) == len(names), key
        assert names, key


@pytest.mark.parametrize("server_kind", ["exported", "anndata"])
def test_the_capability_route_answers_without_credentials_like_the_other_reads(
    tmp_path: Path,
    server_kind: str,
) -> None:
    with _both_servers(tmp_path) as servers:
        address = servers[server_kind]
        health_status, _body, _type, _length = _request(address, "/_cellucid/health")
        anonymous_status, anonymous_body, _type, _length = _request(address, PROTOCOL_PATH)
        # A token is not part of this route's contract, so offering one must not
        # buy a different answer -- and must not be read as authentication.
        credentialed_status, credentialed_body, _type, _length = _request(
            address,
            PROTOCOL_PATH,
            headers={"Authorization": "Bearer not-a-cellucid-credential"},
        )

    assert anonymous_status == health_status == 200
    assert credentialed_status == 200
    assert json.loads(anonymous_body) == json.loads(credentialed_body)
    # The document is a constant of the installed package: no viewer id, no
    # token, no dataset, nothing that could leak by being served anonymously.
    assert json.loads(anonymous_body) == wire._protocol_capability_document()


@pytest.mark.parametrize("server_kind", ["exported", "anndata"])
def test_the_capability_route_describes_its_body_for_head_without_sending_it(
    tmp_path: Path,
    server_kind: str,
) -> None:
    with _both_servers(tmp_path) as servers:
        address = servers[server_kind]
        _status, get_body, _type, _length = _request(address, PROTOCOL_PATH)
        head_status, head_body, head_type, head_length = _request(
            address,
            PROTOCOL_PATH,
            method="HEAD",
        )

    assert head_status == 200
    assert head_body == b""
    assert head_type == "application/json"
    assert head_length == str(len(get_body))


@pytest.mark.parametrize("server_kind", ["exported", "anndata"])
@pytest.mark.parametrize(
    ("method", "path"),
    [
        ("POST", PROTOCOL_PATH),
        ("GET", "/_cellucid/protocol/"),
        ("GET", "/_cellucid/protocols"),
        ("GET", "/_cellucid/Protocol"),
        ("GET", "/_cellucid/protocol/../datasets"),
        ("GET", "/protocol"),
        ("HEAD", "/_cellucid/protocol/"),
    ],
)
def test_only_the_exact_capability_route_is_answered(
    tmp_path: Path,
    server_kind: str,
    method: str,
    path: str,
) -> None:
    with _both_servers(tmp_path) as servers:
        status, body, _type, _length = _request(
            servers[server_kind],
            path,
            method=method,
        )

    assert status == 404, (server_kind, method, path, status, body)
    assert b"events" not in body
