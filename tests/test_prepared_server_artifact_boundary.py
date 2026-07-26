from __future__ import annotations

import http.client
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data import prepare
from cellucid.server import CellucidServer


def _prepare_fixture(output: Path, *, compression: int | None = None) -> None:
    embedding = np.array(
        [[0.0, 0.0], [1.0, 0.25], [0.25, 1.0]],
        dtype=np.float32,
    )
    prepare(
        latent_space=embedding.copy(),
        obs=pd.DataFrame({"score": [0.25, 0.5, 0.75]}),
        X_umap_2d=embedding,
        out_dir=output,
        dataset_name="Server artifact boundary",
        dataset_id="server-artifact-boundary",
        obs_categorical_dtype="uint16",
        centroid_min_points=1,
        compression=compression,
    )


def _request(
    url: str,
    *,
    method: str = "GET",
    headers: dict[str, str] | None = None,
) -> tuple[int, bytes, urllib.response.addinfourl]:
    request = urllib.request.Request(
        url,
        method=method,
        headers={} if headers is None else headers,
    )
    try:
        response = urllib.request.urlopen(request, timeout=5)
    except urllib.error.HTTPError as error:
        return error.code, error.read(), error
    with response:
        return response.status, response.read(), response


def test_prepared_server_serves_only_declared_in_root_regular_artifacts(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    _prepare_fixture(output)
    secret = output / "private-notes.txt"
    secret.write_text("must never leave the server", encoding="utf-8")
    outside = tmp_path / "outside-secret.txt"
    outside.write_text("outside secret", encoding="utf-8")
    symlink = output / "linked-secret.txt"
    symlink.symlink_to(outside)

    server = CellucidServer(
        output,
        port=0,
        open_browser=False,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    try:
        identity_bytes = (output / "dataset_identity.json").read_bytes()
        status, body, headers = _request(f"{server.url}/dataset_identity.json")
        assert status == 200
        assert body == identity_bytes
        assert headers.headers["Content-Type"] == "application/json"

        status, body, headers = _request(
            f"{server.url}/dataset_identity.json",
            method="HEAD",
        )
        assert status == 200
        assert body == b""
        assert int(headers.headers["Content-Length"]) == len(identity_bytes)

        for path in (
            "private-notes.txt",
            "linked-secret.txt",
            "obs/",
            "dataset_identity.json?download=1",
            "%2e%2e/outside-secret.txt",
            "%252e%252e/outside-secret.txt",
        ):
            for method in ("GET", "HEAD"):
                status, body, _headers = _request(
                    f"{server.url}/{path}",
                    method=method,
                )
                assert status == 404, (method, path, status, body)

        point_bytes = (output / "points_2d.bin").read_bytes()
        status, body, headers = _request(
            f"{server.url}/points_2d.bin",
            headers={"Range": "bytes=4-11"},
        )
        assert status == 206
        assert body == point_bytes[4:12]
        assert headers.headers["Content-Range"] == f"bytes 4-11/{len(point_bytes)}"
        assert headers.headers["Accept-Ranges"] == "bytes"

        status, body, headers = _request(
            f"{server.url}/points_2d.bin",
            method="HEAD",
            headers={"Range": "bytes=4-11"},
        )
        assert status == 206
        assert body == b""
        assert int(headers.headers["Content-Length"]) == 8

        status, _body, headers = _request(
            f"{server.url}/points_2d.bin",
            headers={"Range": "bytes=999-1000"},
        )
        assert status == 416
        assert headers.headers["Content-Range"] == f"bytes */{len(point_bytes)}"

        status, body, _headers = _request(f"{server.url}/_cellucid/health")
        assert status == 200
        assert b'"status": "ok"' in body
        status, body, _headers = _request(
            f"{server.url}/_cellucid/health",
            method="HEAD",
        )
        assert status == 200
        assert body == b""
    finally:
        server.stop()


@pytest.mark.parametrize(
    "range_value",
    [
        "bytes=",
        "items=0-1",
        "bytes=0-1,2-3",
        "bytes=9-4",
        "bytes=-0",
        "bytes=999-",
        "bytes=000000000000000000000-1",
    ],
)
def test_prepared_server_rejects_every_noncanonical_or_unsatisfied_range(
    tmp_path: Path,
    range_value: str,
) -> None:
    output = tmp_path / "generation"
    _prepare_fixture(output)
    point_bytes = (output / "points_2d.bin").read_bytes()
    server = CellucidServer(
        output,
        port=0,
        open_browser=False,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    try:
        for method in ("GET", "HEAD"):
            status, body, headers = _request(
                f"{server.url}/points_2d.bin",
                method=method,
                headers={"Range": range_value},
            )
            assert status == 416
            assert headers.headers["Content-Range"] == f"bytes */{len(point_bytes)}"
            assert body == (
                b"Requested byte range is not satisfiable"
                if method == "GET"
                else b""
            )
    finally:
        server.stop()


@pytest.mark.parametrize(
    ("range_value", "expected"),
    [
        ("bytes=4-", slice(4, None)),
        ("bytes=-4", slice(-4, None)),
        ("bytes=4-999", slice(4, None)),
        ("bytes=0-0", slice(0, 1)),
    ],
)
def test_prepared_server_supports_one_exact_byte_range(
    tmp_path: Path,
    range_value: str,
    expected: slice,
) -> None:
    output = tmp_path / "generation"
    _prepare_fixture(output)
    point_bytes = (output / "points_2d.bin").read_bytes()
    server = CellucidServer(
        output,
        port=0,
        open_browser=False,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    try:
        status, body, headers = _request(
            f"{server.url}/points_2d.bin",
            headers={"Range": range_value},
        )
        assert status == 206
        assert body == point_bytes[expected]
        assert int(headers.headers["Content-Length"]) == len(body)
        assert headers.headers["Accept-Ranges"] == "bytes"
    finally:
        server.stop()


def test_prepared_server_rejects_an_artifact_changed_after_validation(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    _prepare_fixture(output)
    server = CellucidServer(
        output,
        port=0,
        open_browser=False,
        quiet=True,
        serve_web_ui=False,
    )
    identity = output / "dataset_identity.json"
    identity.write_bytes(identity.read_bytes() + b"\n")
    server.start_background()
    try:
        for method in ("GET", "HEAD"):
            status, body, _headers = _request(
                f"{server.url}/dataset_identity.json",
                method=method,
            )
            assert status == 409
            assert body == (
                b"Prepared artifact changed after server validation"
                if method == "GET"
                else b""
            )
    finally:
        server.stop()


def test_prepared_server_sends_gzip_artifacts_as_opaque_declared_bytes(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    _prepare_fixture(output, compression=1)
    compressed_path = output / "points_2d.bin.gz"
    compressed_bytes = compressed_path.read_bytes()
    server = CellucidServer(
        output,
        port=0,
        open_browser=False,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    try:
        status, body, headers = _request(f"{server.url}/points_2d.bin.gz")
        assert status == 200
        assert body == compressed_bytes
        assert headers.headers.get("Content-Encoding") is None
        assert int(headers.headers["Content-Length"]) == len(compressed_bytes)
    finally:
        server.stop()


def test_prepared_server_rejects_duplicate_range_headers(tmp_path: Path) -> None:
    output = tmp_path / "generation"
    _prepare_fixture(output)
    point_bytes = (output / "points_2d.bin").read_bytes()
    server = CellucidServer(
        output,
        port=0,
        open_browser=False,
        quiet=True,
        serve_web_ui=False,
    )
    server.start_background()
    try:
        connection = http.client.HTTPConnection("127.0.0.1", server.port, timeout=5)
        try:
            connection.putrequest("GET", "/points_2d.bin")
            connection.putheader("Range", "bytes=0-1")
            connection.putheader("Range", "bytes=2-3")
            connection.endheaders()
            response = connection.getresponse()
            assert response.status == 416
            assert response.getheader("Content-Range") == f"bytes */{len(point_bytes)}"
            assert response.read() == b"Requested byte range is not satisfiable"
        finally:
            connection.close()
    finally:
        server.stop()


def test_prepared_server_rejects_a_declared_symlink_before_binding(
    tmp_path: Path,
) -> None:
    output = tmp_path / "generation"
    _prepare_fixture(output)
    point_path = output / "points_2d.bin"
    outside = tmp_path / "outside.bin"
    outside.write_bytes(point_path.read_bytes())
    point_path.unlink()
    point_path.symlink_to(outside)

    with pytest.raises(ValueError, match="symbolic link"):
        CellucidServer(
            output,
            port=0,
            open_browser=False,
            quiet=True,
            serve_web_ui=False,
        )
