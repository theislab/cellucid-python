"""The viewer server publishes its channel and the capabilities it accepts.

``GET /_cellucid/protocol`` exists so a web build can find the viewer channel
and discover whether the notebook driving it accepts an event before it emits
one. The route therefore
has to answer without credentials
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
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import cellucid.jupyter._wire as wire
from cellucid.anndata_server import AnnDataServer

PROTOCOL_PATH = "/_cellucid/protocol"


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


def test_the_server_publishes_the_declaration_the_validator_reads() -> None:
    expected = wire._protocol_capability_document()
    with _anndata_server() as address:
        status, body, content_type, _length = _request(
            address,
            PROTOCOL_PATH,
        )

    assert status == 200
    assert content_type == "application/json"
    # Compared against the module the validator itself reads, never a literal:
    # a copy kept in the server is exactly the drift this route must not have.
    assert json.loads(body) == expected


def test_the_capability_document_carries_only_growable_capability_lists() -> None:
    with _anndata_server() as address:
        _status, body, _type, _length = _request(address, PROTOCOL_PATH)
    document: Any = json.loads(body)

    # The key set is fixed on purpose. A viewer validates payloads against exact
    # key sets, so a later capability has to arrive inside these lists rather
    # than as a new key, which would break the builds this route exists to
    # protect. There is deliberately no version to compare against either.
    # `binaryCommands` is the one key added since: commands that name cells
    # became binary frames, a form an older viewer cannot read, so an older
    # viewer refuses this host at once rather than misread its commands.
    assert set(document) == {"channel", "commands", "binaryCommands", "events"}
    assert set(document["binaryCommands"]) <= set(document["commands"])
    # Where the viewer opens its channel, relative to the server's address.
    assert document.pop("channel") == "_cellucid/ws"
    for key, names in document.items():
        assert isinstance(names, list), key
        assert all(type(name) is str for name in names), key
        assert names == sorted(names), key
        assert len(set(names)) == len(names), key
        assert names, key


def test_the_capability_route_answers_without_credentials_like_the_other_reads() -> None:
    with _anndata_server() as address:
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


def test_the_capability_route_describes_its_body_for_head_without_sending_it() -> None:
    with _anndata_server() as address:
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
    method: str,
    path: str,
) -> None:
    with _anndata_server() as address:
        status, body, _type, _length = _request(
            address,
            path,
            method=method,
        )

    # This server takes no POST at all; every other near miss is not a route.
    assert status == (405 if method == "POST" else 404), (method, path, status, body)
    assert b"events" not in body
