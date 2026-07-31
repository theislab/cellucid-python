"""Host-header contract: the servers answer only their own loopback authority.

A page on any origin can publish a short-lived DNS record that re-resolves its
own hostname to ``127.0.0.1``. The browser then treats the local Cellucid data
server as same-origin, so no CORS check applies and the response body is handed
to the attacker's page. The bind address stops remote clients but not rebinding;
only ``Host`` validation does.
"""

from __future__ import annotations

import json
import logging
import socket
import urllib.request
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from http.client import HTTPMessage
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellucid._server_base import (
    CORSMixin,
    _host_authority_is_loopback,
    _parse_host_authority,
    _requires_loopback_host_header,
    require_allowed_hosts,
)
from cellucid.anndata_server import AnnDataServer
from cellucid.server import CellucidServer

REBOUND_HOST = "cellucid-rebind.evil.example"
# The authority a jupyter-server-proxy front end forwards verbatim.
PROXY_HOST = "hub.example.org"
POINTS_BYTES = b"\x11\x22\x33\x44\x55\x66\x77\x88"


def _prepared_dataset(root: Path) -> Path:
    root.mkdir()
    (root / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": "host-contract",
                "name": "Host contract",
                "description": "",
                "stats": {
                    "n_cells": 1,
                    "n_genes": 0,
                    "n_obs_fields": 0,
                    "n_categorical_fields": 0,
                    "n_continuous_fields": 0,
                    "has_connectivity": False,
                    "n_edges": None,
                },
                "embeddings": {
                    "available_dimensions": [2],
                    "default_dimension": 2,
                    "files": {"2d": "points_2d.bin"},
                },
                "obs_fields": [],
                "export_settings": {
                    "compression": None,
                    "var_quantization": None,
                    "obs_continuous_quantization": None,
                    "obs_categorical_dtype": "uint16",
                },
            }
        ),
        encoding="utf-8",
    )
    (root / "obs_manifest.json").write_text(
        json.dumps(
            {
                "_format": "compact_v1",
                "n_points": 1,
                "centroid_outlier_quantile": None,
                "latent_key": "latent_space",
                "compression": None,
                "_obsSchemas": {},
                "_continuousFields": [],
                "_categoricalFields": [],
            }
        ),
        encoding="utf-8",
    )
    (root / "points_2d.bin").write_bytes(POINTS_BYTES)
    return root


def _anndata() -> AnnData:
    return AnnData(
        X=np.array([[1.0], [2.0]], dtype=np.float32),
        obs=pd.DataFrame(
            {"label": pd.Categorical(["A", "B"])},
            index=["cell-1", "cell-2"],
        ),
        var=pd.DataFrame(index=["gene-1"]),
        obsm={"X_umap_2d": np.array([[0.0, 0.0], [1.0, 1.0]], dtype=np.float32)},
    )


@contextmanager
def _prepared_server(
    tmp_path: Path,
    *,
    allowed_hosts: Sequence[str] | None = None,
) -> Iterator[tuple[CellucidServer, int]]:
    server = CellucidServer(
        _prepared_dataset(tmp_path / "dataset"),
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        allowed_hosts=allowed_hosts,
    )
    server.start_background()
    try:
        yield server, server.port
    finally:
        server.stop()


@contextmanager
def _anndata_server(
    *,
    allowed_hosts: Sequence[str] | None = None,
) -> Iterator[tuple[AnnDataServer, int]]:
    server = AnnDataServer(
        _anndata(),
        dataset_name="Host contract",
        dataset_id="host-contract",
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=False,
        allowed_hosts=allowed_hosts,
    )
    server.start_background()
    try:
        yield server, server.port
    finally:
        server.stop()


def _exchange(port: int, request: bytes) -> tuple[int, bytes]:
    """Send one exact request over a raw socket and read the whole response."""
    with socket.create_connection(("127.0.0.1", port), timeout=5) as connection:
        connection.sendall(request)
        chunks: list[bytes] = []
        while True:
            chunk = connection.recv(65536)
            if not chunk:
                break
            chunks.append(chunk)
    raw = b"".join(chunks)
    status_line = raw.split(b"\r\n", 1)[0]
    fields = status_line.split(b" ")
    assert fields[0].startswith(b"HTTP/"), raw
    return int(fields[1]), raw


def _request(
    target: str,
    host_values: Sequence[str],
    *,
    method: str = "GET",
) -> bytes:
    lines = [f"{method} {target} HTTP/1.1"]
    lines.extend(f"Host: {value}" for value in host_values)
    lines.append("Connection: close")
    return ("\r\n".join(lines) + "\r\n\r\n").encode("ascii")


# ---------------------------------------------------------------------------
# The attack
# ---------------------------------------------------------------------------


def test_rebound_dns_name_never_reaches_a_route_on_the_prepared_server(
    tmp_path: Path,
) -> None:
    """A rebound hostname must not be served prepared artifacts."""
    with _prepared_server(tmp_path) as (_server, port):
        status, raw = _exchange(
            port,
            _request("/points_2d.bin", [f"{REBOUND_HOST}:{port}"]),
        )

    assert status == 421
    # The private payload was never routed, let alone written.
    assert POINTS_BYTES not in raw
    # The attacker-supplied authority is never reflected back.
    assert REBOUND_HOST.encode("ascii") not in raw


def test_rebound_dns_name_never_reaches_a_route_on_the_anndata_server() -> None:
    """One shared check covers the AnnData server too."""
    with _anndata_server() as (_server, port):
        status, raw = _exchange(
            port,
            _request("/dataset_identity.json", [f"{REBOUND_HOST}:{port}"]),
        )

    assert status == 421
    assert b"host-contract" not in raw
    assert REBOUND_HOST.encode("ascii") not in raw


def test_rebound_dns_name_is_refused_on_every_method(tmp_path: Path) -> None:
    """The check precedes dispatch, so no verb can slip past it."""
    with _prepared_server(tmp_path) as (_server, port):
        results = {
            method: _exchange(
                port,
                _request(
                    "/_cellucid/events" if method == "POST" else "/_cellucid/info",
                    [REBOUND_HOST],
                    method=method,
                ),
            )[0]
            for method in ("GET", "HEAD", "POST", "OPTIONS")
        }

    assert results == {"GET": 421, "HEAD": 421, "POST": 421, "OPTIONS": 421}


@pytest.mark.parametrize(
    "host_values",
    [
        pytest.param([], id="missing"),
        pytest.param([""], id="empty"),
        pytest.param(["127.0.0.1", "127.0.0.1"], id="duplicated"),
        pytest.param(["127.0.0.1"], id="loopback-with-nondefault-port-elided"),
        pytest.param(["127.0.0.1:1"], id="wrong-port"),
        pytest.param(["localhost:1"], id="wrong-port-by-name"),
        pytest.param(["::1"], id="unbracketed-ipv6"),
        pytest.param(["evil@127.0.0.1"], id="userinfo"),
        pytest.param(["127.0.0.1/evil"], id="path"),
        pytest.param(["127.0.0.1:"], id="empty-port"),
        pytest.param(["localhost.evil.example"], id="localhost-prefixed-name"),
        pytest.param(["localhost.local"], id="localhost-suffixed-name"),
        pytest.param(["0.0.0.0"], id="unspecified-address"),
        pytest.param(["10.0.0.1"], id="private-address"),
    ],
)
def test_malformed_or_foreign_host_authorities_are_refused(
    tmp_path: Path,
    host_values: list[str],
) -> None:
    with _prepared_server(tmp_path) as (_server, port):
        status, _raw = _exchange(port, _request("/_cellucid/info", host_values))

    assert status == 421


def test_a_refusal_is_reported_to_the_operator_and_not_to_the_client(
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Requests log at DEBUG, so a refusal would otherwise be undiagnosable."""
    with (
        caplog.at_level(logging.WARNING, logger="cellucid.server"),
        _prepared_server(tmp_path) as (_server, port),
    ):
        status, raw = _exchange(
            port,
            _request("/_cellucid/info", [f"{REBOUND_HOST}:{port}"]),
        )

    assert status == 421
    assert REBOUND_HOST.encode("ascii") not in raw
    warnings = [
        record.getMessage() for record in caplog.records if record.levelno >= logging.WARNING
    ]
    assert any(REBOUND_HOST in message for message in warnings)


def test_a_duplicated_correct_host_is_still_refused(tmp_path: Path) -> None:
    """Two Host headers are ambiguous to intermediaries and never legitimate."""
    with _prepared_server(tmp_path) as (_server, port):
        authority = f"127.0.0.1:{port}"
        status, _raw = _exchange(port, _request("/_cellucid/info", [authority, authority]))

    assert status == 421


# ---------------------------------------------------------------------------
# Legitimate local use
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "template",
    ["127.0.0.1:{port}", "localhost:{port}", "[::1]:{port}", "LOCALHOST:{port}"],
)
def test_loopback_authorities_are_served(tmp_path: Path, template: str) -> None:
    with _prepared_server(tmp_path) as (_server, port):
        status, raw = _exchange(
            port,
            _request("/points_2d.bin", [template.format(port=port)]),
        )

    assert status == 200
    assert raw.endswith(POINTS_BYTES)


def test_the_real_serving_path_still_works_end_to_end(tmp_path: Path) -> None:
    """The normal client path — the one the app and `jupyter.py` both use."""
    with _prepared_server(tmp_path) as (server, _port):
        base = server.url
        with urllib.request.urlopen(f"{base}/_cellucid/health", timeout=5) as response:
            health = json.loads(response.read())
        with urllib.request.urlopen(f"{base}/dataset_identity.json", timeout=5) as response:
            identity = json.loads(response.read())
        with urllib.request.urlopen(f"{base}/points_2d.bin", timeout=5) as response:
            points = response.read()

    assert health["status"] == "ok"
    assert identity["id"] == "host-contract"
    assert points == POINTS_BYTES


def test_the_anndata_serving_path_still_works_end_to_end() -> None:
    with _anndata_server() as (server, _port):
        base = server.url
        with urllib.request.urlopen(f"{base}/_cellucid/health", timeout=5) as response:
            health = json.loads(response.read())

    assert health["status"] == "ok"
    assert health["type"] == "anndata"


# ---------------------------------------------------------------------------
# Authority parsing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("127.0.0.1:8765", ("127.0.0.1", 8765)),
        ("localhost:8765", ("localhost", 8765)),
        ("localhost", ("localhost", None)),
        ("[::1]:8765", ("::1", 8765)),
        ("[::1]", ("::1", None)),
        ("[0:0:0:0:0:0:0:1]:8765", ("0:0:0:0:0:0:0:1", 8765)),
    ],
)
def test_host_authorities_are_parsed_without_naive_colon_splitting(
    value: str,
    expected: tuple[str, int | None],
) -> None:
    assert _parse_host_authority(value) == expected


@pytest.mark.parametrize(
    "value",
    [
        "",
        "::1",
        "[::1",
        "[::1]x",
        "[::1]:",
        "[fe80::1%eth0]:8765",
        "[not-an-address]:8765",
        "127.0.0.1:08765",
        "127.0.0.1:0",
        "127.0.0.1:65536",
        "127.0.0.1:87 65",
        "127.0.0.1:87.65",
        "127.0.0.1:-1",
        "user@127.0.0.1:8765",
        "127.0.0.1:8765/path",
        ":8765",
        "127.0.0.1:8765:8765",
        "127.0.0.1\t",
        "localhost\n",
        "lócalhost",
        b"127.0.0.1:8765",
        None,
    ],
)
def test_malformed_host_authorities_are_rejected_by_the_parser(value: object) -> None:
    assert _parse_host_authority(value) is None


@pytest.mark.parametrize(
    ("host", "expected"),
    [
        ("127.0.0.1", True),
        ("127.0.0.2", True),
        ("127.255.255.254", True),
        ("::1", True),
        ("0:0:0:0:0:0:0:1", True),
        ("localhost", True),
        ("LocalHost", True),
        ("0.0.0.0", False),
        ("10.0.0.1", False),
        ("192.168.1.5", False),
        ("::", False),
        ("localhost.evil.example", False),
        ("evil.example", False),
        # RFC 6761 reserves the name, not names that merely contain it.
        ("localhost.", False),
        ("cellucid.com", False),
    ],
)
def test_only_loopback_literals_and_localhost_are_authorized(host: str, expected: bool) -> None:
    assert _host_authority_is_loopback(host) is expected


@pytest.mark.parametrize(
    ("bound_host", "expected"),
    [
        # Loopback binds — the default — are enforced.
        ("127.0.0.1", True),
        ("127.0.0.5", True),
        ("::1", True),
        # A deliberate non-loopback bind is an explicit request for network
        # exposure, and the legitimate names of an exposed deployment cannot be
        # known here, so the loopback rule is not imposed on it.
        ("0.0.0.0", False),
        ("::", False),
        ("192.168.1.5", False),
        # Anything undetermined is enforced rather than waved through.
        ("not-an-address", True),
        (b"127.0.0.1", True),
        (None, True),
    ],
)
def test_enforcement_follows_the_actual_bound_address(bound_host: object, expected: bool) -> None:
    assert _requires_loopback_host_header(bound_host) is expected


def test_enforcement_uses_the_os_assigned_port_not_the_requested_one(
    tmp_path: Path,
) -> None:
    """Port 0 is a request for an assigned port; the check must follow it."""
    with _prepared_server(tmp_path) as (server, port):
        assert server._server is not None
        assert server._server.server_port == port
        assert port != 0
        accepted, _raw = _exchange(port, _request("/_cellucid/info", [f"127.0.0.1:{port}"]))
        refused, _raw = _exchange(port, _request("/_cellucid/info", ["127.0.0.1:8765"]))

    assert accepted == 200
    assert refused == 421


# ---------------------------------------------------------------------------
# The declared reverse-proxy authority (`allowed_hosts`)
# ---------------------------------------------------------------------------
#
# `jupyter-server-proxy` copies the incoming request headers verbatim, so a
# loopback-bound Cellucid server behind it receives the hub's own `Host`. That
# is byte-for-byte what a rebound page sends, and a rebound page is same-origin
# so it can forge any `X-Forwarded-*` header it likes: no in-band signal can
# tell the two apart. The proxy can only be declared, never detected.


def test_a_declared_proxy_authority_is_served_the_real_payload(tmp_path: Path) -> None:
    """The documented `jupyter-server-proxy` path works again, end to end."""
    with _prepared_server(tmp_path, allowed_hosts=[PROXY_HOST]) as (_server, port):
        status, raw = _exchange(port, _request("/points_2d.bin", [PROXY_HOST]))

    assert status == 200
    assert raw.endswith(POINTS_BYTES)


def test_a_declared_proxy_authority_is_served_by_the_anndata_server() -> None:
    with _anndata_server(allowed_hosts=[PROXY_HOST]) as (_server, port):
        status, raw = _exchange(port, _request("/dataset_identity.json", [PROXY_HOST]))

    assert status == 200
    assert b"host-contract" in raw


@pytest.mark.parametrize(
    "template",
    [
        # A proxy publishes its own front-end port, which has nothing to do with
        # the loopback port bound here, so a declared name matches on any port.
        "hub.example.org",
        "hub.example.org:80",
        "hub.example.org:443",
        "hub.example.org:8000",
        # Authority hosts are case-insensitive.
        "HUB.Example.ORG:8000",
        # ...including the port this server happens to be bound on.
        "hub.example.org:{port}",
    ],
)
def test_a_declared_name_matches_on_whatever_port_the_proxy_publishes(
    tmp_path: Path,
    template: str,
) -> None:
    with _prepared_server(tmp_path, allowed_hosts=[PROXY_HOST]) as (_server, port):
        status, _raw = _exchange(
            port,
            _request("/_cellucid/info", [template.format(port=port)]),
        )

    assert status == 200


def test_declaring_one_name_does_not_admit_any_other(tmp_path: Path) -> None:
    """The hole stays shut for every authority the operator did not name."""
    with _prepared_server(tmp_path, allowed_hosts=[PROXY_HOST]) as (_server, port):
        results = {
            value: _exchange(port, _request("/points_2d.bin", [value]))
            for value in (
                f"{REBOUND_HOST}:{port}",
                REBOUND_HOST,
                # A prefix/suffix of a declared name is a different name.
                f"evil.{PROXY_HOST}",
                f"{PROXY_HOST}.evil.example",
                "hub.example.orgx",
                # A wildcard is never expanded, so it can never be a match.
                "*.example.org",
            )
        }

    assert {value: status for value, (status, _raw) in results.items()} == {
        f"{REBOUND_HOST}:{port}": 421,
        REBOUND_HOST: 421,
        f"evil.{PROXY_HOST}": 421,
        f"{PROXY_HOST}.evil.example": 421,
        "hub.example.orgx": 421,
        "*.example.org": 421,
    }
    for _status, raw in results.values():
        assert POINTS_BYTES not in raw


def test_the_proxy_authority_is_refused_without_the_declaration(tmp_path: Path) -> None:
    """No opt-in means byte-identical behaviour to the loopback-only guard."""
    with _prepared_server(tmp_path) as (_server, port):
        bare, _raw = _exchange(port, _request("/points_2d.bin", [PROXY_HOST]))
        ported, _raw = _exchange(port, _request("/points_2d.bin", [f"{PROXY_HOST}:{port}"]))

    assert bare == 421
    assert ported == 421


def test_the_proxy_authority_is_refused_without_the_declaration_on_anndata() -> None:
    with _anndata_server() as (_server, port):
        status, raw = _exchange(port, _request("/dataset_identity.json", [PROXY_HOST]))

    assert status == 421
    assert b"host-contract" not in raw


def test_the_local_loopback_authority_still_works_beside_a_declared_name(
    tmp_path: Path,
) -> None:
    """Declaring a proxy name must not cost the operator their own `curl`."""
    with _prepared_server(tmp_path, allowed_hosts=[PROXY_HOST]) as (server, port):
        with urllib.request.urlopen(f"{server.url}/points_2d.bin", timeout=5) as response:
            points = response.read()
        wrong_port, _raw = _exchange(port, _request("/_cellucid/info", ["127.0.0.1:8765"]))

    assert points == POINTS_BYTES
    # The loopback branch still pins the bound port; only declared names do not.
    assert wrong_port == 421


def test_a_declaration_does_not_relax_the_single_host_header_rule(tmp_path: Path) -> None:
    with _prepared_server(tmp_path, allowed_hosts=[PROXY_HOST]) as (_server, port):
        duplicated, _raw = _exchange(port, _request("/_cellucid/info", [PROXY_HOST, PROXY_HOST]))
        malformed, _raw = _exchange(port, _request("/_cellucid/info", [f"user@{PROXY_HOST}"]))
        missing, _raw = _exchange(port, _request("/_cellucid/info", []))

    assert duplicated == 421
    assert malformed == 421
    assert missing == 421


def test_a_declared_loopback_literal_is_matched_in_canonical_form(tmp_path: Path) -> None:
    """An entry and a header are reduced through one normalization, not two."""
    with _prepared_server(tmp_path, allowed_hosts=["::1"]) as (_server, port):
        expanded, _raw = _exchange(
            port,
            _request("/_cellucid/info", ["[0:0:0:0:0:0:0:1]:1"]),
        )

    assert expanded == 200


# ---------------------------------------------------------------------------
# Where a declaration applies (bind address)
# ---------------------------------------------------------------------------


def _authorization(
    *,
    bound_host: object,
    bound_port: int,
    host_values: Sequence[str],
    allowed_hosts: Sequence[str] | None,
) -> bool:
    """Drive the exact Host decision without binding a public interface."""

    class _Probe(CORSMixin):
        def _request_server(self):  # type: ignore[override]
            return SimpleNamespace(
                server_address=(bound_host, bound_port),
                server_port=bound_port,
            )

        def _request_headers(self):  # type: ignore[override]
            message = HTTPMessage()
            for value in host_values:
                message["Host"] = value
            return message

    probe = _Probe()
    probe.allowed_hosts = require_allowed_hosts(allowed_hosts)
    return probe._host_header_is_authorized()


@pytest.mark.parametrize(
    ("bound_host", "host_value", "allowed_hosts", "expected"),
    [
        # Unchanged: a non-loopback bind with no declaration is not enforced.
        ("0.0.0.0", REBOUND_HOST, None, True),
        ("::", REBOUND_HOST, None, True),
        ("192.168.1.5", REBOUND_HOST, None, True),
        # A declaration is the operator stating the complete set of legitimate
        # names, which is exactly the knowledge a non-loopback bind lacked. It
        # can only be supplied deliberately, so honouring it there is a strict
        # tightening no existing deployment can be surprised by.
        ("0.0.0.0", PROXY_HOST, [PROXY_HOST], True),
        ("0.0.0.0", REBOUND_HOST, [PROXY_HOST], False),
        ("192.168.1.5", REBOUND_HOST, [PROXY_HOST], False),
        # The operator's own loopback access survives an exposed bind too.
        ("0.0.0.0", "127.0.0.1:8765", [PROXY_HOST], True),
        ("0.0.0.0", "127.0.0.1:1", [PROXY_HOST], False),
        # Loopback binds behave the same with and without a declaration.
        ("127.0.0.1", "localhost:8765", None, True),
        ("127.0.0.1", "localhost:8765", [PROXY_HOST], True),
        ("127.0.0.1", REBOUND_HOST, None, False),
        ("127.0.0.1", REBOUND_HOST, [PROXY_HOST], False),
    ],
)
def test_a_declaration_governs_where_the_check_applies(
    bound_host: str,
    host_value: str,
    allowed_hosts: list[str] | None,
    expected: bool,
) -> None:
    assert (
        _authorization(
            bound_host=bound_host,
            bound_port=8765,
            host_values=[host_value],
            allowed_hosts=allowed_hosts,
        )
        is expected
    )


def test_a_handler_without_a_validated_declaration_fails_closed() -> None:
    """A missing or wrong-typed attribute must mean 'strictest', never 'off'."""

    class _Probe(CORSMixin):
        def _request_server(self):  # type: ignore[override]
            return SimpleNamespace(server_address=("127.0.0.1", 8765), server_port=8765)

        def _request_headers(self):  # type: ignore[override]
            message = HTTPMessage()
            message["Host"] = PROXY_HOST
            return message

    unset = _Probe()
    assert unset._allowed_host_names() == frozenset()
    assert unset._host_header_is_authorized() is False

    forged = _Probe()
    # A plain set, a list or a string is not what `require_allowed_hosts`
    # returns, so it is not a declaration.
    for value in ({PROXY_HOST}, [PROXY_HOST], PROXY_HOST, None):
        forged.allowed_hosts = value  # type: ignore[assignment]
        assert forged._allowed_host_names() == frozenset()
        assert forged._host_header_is_authorized() is False


# ---------------------------------------------------------------------------
# Validating the declaration itself
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (None, frozenset()),
        ([], frozenset()),
        ((), frozenset()),
        ([PROXY_HOST], frozenset({PROXY_HOST})),
        # Names are case-insensitive and de-duplicated.
        ([PROXY_HOST, "HUB.Example.ORG"], frozenset({PROXY_HOST})),
        # IP literals are canonicalized, so one spelling cannot hide another.
        (["0:0:0:0:0:0:0:1"], frozenset({"::1"})),
        (["::1", "0:0:0:0:0:0:0:1"], frozenset({"::1"})),
        (["127.0.0.1"], frozenset({"127.0.0.1"})),
        (["a-b.c-d.example", "x1.example"], frozenset({"a-b.c-d.example", "x1.example"})),
    ],
)
def test_a_declaration_is_reduced_to_exact_comparison_form(
    value: object,
    expected: frozenset[str],
) -> None:
    assert require_allowed_hosts(value) == expected


@pytest.mark.parametrize(
    ("entry", "message"),
    [
        ("", "is empty"),
        (" ", "contains whitespace"),
        ("hub.example.org ", "contains whitespace"),
        ("a b.example", "contains whitespace"),
        ("hub.example.org\n", "contains whitespace"),
        ("*", "wildcard"),
        ("*.example.org", "wildcard"),
        ("hub.*.org", "wildcard"),
        ("http://hub.example.org", "URL"),
        ("hub.example.org/path", "URL"),
        ("user@hub.example.org", "credentials"),
        ("[::1]", "bracketed"),
        ("[::1]:8000", "bracketed"),
        ("fe80::1%eth0", "zone identifier"),
        ("hub.example.org:8000", "carries a port"),
        ("127.0.0.1:8765", "carries a port"),
        ("hub.example.org.", "ends with a dot"),
        ("localhost.", "ends with a dot"),
        ("-hub.example.org", "is not a host name"),
        ("hub-.example.org", "is not a host name"),
        ("hub..example.org", "is not a host name"),
        (".hub.example.org", "is not a host name"),
        ("hub_example.org", "is not a host name"),
        ("a" * 64 + ".example", "is not a host name"),
        ("a" * 250 + ".example", "is not a host name"),
        ("hüb.example.org", "not ASCII"),
        ("hub.example.org\x00", "control character"),
        ("hub.example.org\x7f", "control character"),
    ],
)
def test_every_unsupported_declaration_entry_is_refused_by_name(
    entry: str,
    message: str,
) -> None:
    """A silently dropped or silently widened entry is the hole with a UI."""
    with pytest.raises(ValueError, match=message):
        require_allowed_hosts([entry])


@pytest.mark.parametrize(
    ("value", "message"),
    [
        (PROXY_HOST, "not one string"),
        (b"hub.example.org", "not one string"),
        (5, "list or tuple"),
        ({PROXY_HOST: True}, "list or tuple"),
        ({PROXY_HOST}, "list or tuple"),
        ([5], "entries must be strings"),
        ([None], "entries must be strings"),
        ([b"hub.example.org"], "entries must be strings"),
    ],
)
def test_a_declaration_of_the_wrong_shape_is_refused(value: object, message: str) -> None:
    with pytest.raises(TypeError, match=message):
        require_allowed_hosts(value)


@pytest.mark.parametrize(
    "entry",
    ["*", "hub.example.org:8000", "[::1]", "not a host"],
)
def test_a_bad_declaration_is_refused_at_construction_not_first_request(
    tmp_path: Path,
    entry: str,
) -> None:
    """A typo must fail where the operator can see it, before anything binds."""
    with pytest.raises(ValueError, match="allowed_hosts entry"):
        CellucidServer(
            _prepared_dataset(tmp_path / "dataset"),
            host="127.0.0.1",
            port=0,
            quiet=True,
            serve_web_ui=False,
            allowed_hosts=[entry],
        )

    with pytest.raises(ValueError, match="allowed_hosts entry"):
        AnnDataServer(
            _anndata(),
            dataset_name="Host contract",
            dataset_id="host-contract",
            host="127.0.0.1",
            port=0,
            quiet=True,
            serve_web_ui=False,
            allowed_hosts=[entry],
        )


def test_the_declaration_is_validated_before_the_dataset_is_read(tmp_path: Path) -> None:
    """Ordering proof: a missing directory is never reached for a bad entry."""
    with pytest.raises(ValueError, match="allowed_hosts entry"):
        CellucidServer(
            tmp_path / "does-not-exist",
            host="127.0.0.1",
            port=0,
            quiet=True,
            serve_web_ui=False,
            allowed_hosts=["*"],
        )
