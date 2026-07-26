from __future__ import annotations

import hashlib
import http.client
import io
import json
from pathlib import Path
from unittest import mock

import pytest


def _asset_payloads(build_id: str) -> dict[str, tuple[bytes, str]]:
    return {
        "assets/css/main.css": (
            b"body { color: #123456; }\n",
            "text/css; charset=utf-8",
        ),
        "assets/js/app.js": (
            b"export const ready = true;\n",
            "application/javascript; charset=utf-8",
        ),
        "index.html": (
            (
                "<!doctype html><html><head>"
                f'<meta name="cellucid-web-build-id" content="{build_id}">'
                "</head><body></body></html>"
            ).encode(),
            "text/html; charset=utf-8",
        ),
    }


def _inventory_payload(
    build_id: str,
    assets: dict[str, tuple[bytes, str]],
) -> bytes:
    return json.dumps(
        {
            "version": 1,
            "build_id": build_id,
            "assets": [
                {
                    "path": path,
                    "sha256": hashlib.sha256(payload).hexdigest(),
                    "bytes": len(payload),
                    "content_type": content_type,
                }
                for path, (payload, content_type) in sorted(assets.items())
            ],
        },
        separators=(",", ":"),
    ).encode()


def _install_fetch_fixture(
    monkeypatch,
    *,
    build_id: str,
    assets: dict[str, tuple[bytes, str]],
    corrupt_path: str | None = None,
) -> list[str]:
    from cellucid import web_cache

    inventory = _inventory_payload(build_id, assets)
    calls: list[str] = []

    def fetch(
        *,
        source_url: str,
        asset_path: str,
        timeout: float,
    ):
        assert source_url == "https://viewer.example"
        assert timeout > 0
        calls.append(asset_path)
        if asset_path == web_cache.WEB_ASSET_INVENTORY_FILENAME:
            return web_cache.FetchedWebResponse(
                data=inventory,
                content_type="application/json; charset=utf-8",
                content_length=len(inventory),
            )
        payload, content_type = assets[asset_path]
        if asset_path == corrupt_path:
            payload += b"corrupt"
        return web_cache.FetchedWebResponse(
            data=payload,
            content_type=content_type,
            content_length=len(payload),
        )

    monkeypatch.setattr(web_cache, "_fetch_web_response", fetch)
    return calls


def _cache_snapshot(cache_dir: Path) -> dict[str, bytes]:
    return {
        path.relative_to(cache_dir).as_posix(): path.read_bytes()
        for path in cache_dir.rglob("*")
        if path.is_file()
    }


def _write_current_prepared_dataset(
    root: Path,
    *,
    dataset_id: str,
    dataset_name: str,
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "dataset_identity.json").write_text(
        json.dumps(
            {
                "version": 2,
                "id": dataset_id,
                "name": dataset_name,
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
    (root / "points_2d.bin").write_bytes(b"\x00" * 8)


@pytest.mark.parametrize(
    "payload",
    [
        b"<!doctype html><html><head></head></html>",
        (
            b'<meta name="cellucid-web-build-id" content="one">'
            b'<meta name="cellucid-web-build-id" content="two">'
        ),
        b'<!-- <meta name="cellucid-web-build-id" content="commented"> -->',
    ],
)
def test_index_requires_one_real_build_meta_element(payload: bytes) -> None:
    from cellucid._server_base import _extract_web_build_id

    with pytest.raises(ValueError, match="exactly one"):
        _extract_web_build_id(payload)


@pytest.mark.parametrize(
    "mutate",
    [
        lambda value: {**value, "unexpected": True},
        lambda value: {key: item for key, item in value.items() if key != "version"},
        lambda value: {**value, "version": True},
        lambda value: {**value, "build_id": ""},
        lambda value: {**value, "assets": list(reversed(value["assets"]))},
        lambda value: {
            **value,
            "assets": [*value["assets"], dict(value["assets"][0])],
        },
        lambda value: {
            **value,
            "assets": [
                {**value["assets"][0], "path": "/assets/css/main.css"},
                *value["assets"][1:],
            ],
        },
        lambda value: {
            **value,
            "assets": [
                {**value["assets"][0], "path": "assets/../main.css"},
                *value["assets"][1:],
            ],
        },
        lambda value: {
            **value,
            "assets": [
                {**value["assets"][0], "sha256": "A" * 64},
                *value["assets"][1:],
            ],
        },
        lambda value: {
            **value,
            "assets": [
                {**value["assets"][0], "bytes": -1},
                *value["assets"][1:],
            ],
        },
        lambda value: {
            **value,
            "assets": [
                {**value["assets"][0], "content_type": "TEXT/CSS"},
                *value["assets"][1:],
            ],
        },
        lambda value: {
            **value,
            "assets": [
                asset for asset in value["assets"] if asset["path"] != "index.html"
            ],
        },
        lambda value: {
            **value,
            "assets": [
                {
                    **value["assets"][0],
                    "path": "cellucid-web-assets.json",
                },
                *value["assets"][1:],
            ],
        },
        lambda value: {
            **value,
            "assets": [
                {**value["assets"][0], "path": "assets/css/CON.css"},
                *value["assets"][1:],
            ],
        },
        lambda value: {
            **value,
            "assets": [
                {**value["assets"][0], "path": "assets/css/main:dark.css"},
                *value["assets"][1:],
            ],
        },
    ],
)
def test_inventory_parser_rejects_every_noncanonical_shape(mutate) -> None:
    from cellucid.web_cache import parse_web_asset_inventory

    assets = _asset_payloads("build-1")
    canonical = json.loads(_inventory_payload("build-1", assets))
    invalid = mutate(canonical)

    with pytest.raises((TypeError, ValueError)):
        parse_web_asset_inventory(
            json.dumps(invalid, separators=(",", ":")).encode(),
        )


def test_force_refresh_downloads_one_complete_verified_generation(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    calls = _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
    )

    summary = web_cache.ensure_web_ui_cached(
        cache_dir=tmp_path / "cache",
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )

    assert calls == [
        "cellucid-web-assets.json",
        "assets/css/main.css",
        "assets/js/app.js",
        "index.html",
    ]
    assert summary.build_id == "build-1"
    assert summary.downloaded_files == 3
    assert summary.downloaded_bytes == sum(len(payload) for payload, _ in assets.values())
    assert summary.skipped_files == 0
    assert not hasattr(summary, "errors")
    verified = web_cache.verify_web_ui_cache(
        tmp_path / "cache",
        expected_source_url="https://viewer.example",
    )
    assert verified.build_id == "build-1"
    for path, (payload, _content_type) in assets.items():
        assert (tmp_path / "cache" / path).read_bytes() == payload


def test_prefetch_rejects_internally_inconsistent_inventory_http_metadata(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    inventory = _inventory_payload("build-1", assets)

    monkeypatch.setattr(
        web_cache,
        "_fetch_web_response",
        mock.Mock(
            return_value=web_cache.FetchedWebResponse(
                data=inventory,
                content_type="application/json; charset=utf-8",
                content_length=len(inventory) + 1,
            )
        ),
    )

    with pytest.raises(ValueError, match="body length|Content-Length"):
        web_cache.ensure_web_ui_cached(
            cache_dir=tmp_path / "cache",
            source_url="https://viewer.example",
            force=True,
            show_progress=False,
        )
    assert not (tmp_path / "cache").exists()


def test_complete_source_bound_cache_is_verified_without_network(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
    )
    cache_dir = tmp_path / "cache"
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )

    fetch = mock.Mock(side_effect=AssertionError("cache-only mode used the network"))
    monkeypatch.setattr(
        web_cache,
        "_fetch_web_response",
        fetch,
    )
    summary = web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=False,
        show_progress=False,
    )

    fetch.assert_not_called()
    assert summary.build_id == "build-1"
    assert summary.downloaded_files == 0
    assert summary.downloaded_bytes == 0
    assert summary.skipped_files == 3


@pytest.mark.parametrize("cache_state", ["missing", "empty"])
def test_cache_only_mode_rejects_missing_or_empty_cache_without_network(
    monkeypatch,
    tmp_path: Path,
    cache_state: str,
) -> None:
    from cellucid import web_cache

    cache_dir = tmp_path / "cache"
    if cache_state == "empty":
        cache_dir.mkdir()
    fetch = mock.Mock(side_effect=AssertionError("cache-only mode used the network"))
    monkeypatch.setattr(web_cache, "_fetch_web_response", fetch)

    with pytest.raises(ValueError, match="does not exist|source control is missing"):
        web_cache.ensure_web_ui_cached(
            cache_dir=cache_dir,
            source_url="https://viewer.example",
            force=False,
            show_progress=False,
        )

    fetch.assert_not_called()
    if cache_state == "missing":
        assert not cache_dir.exists()
    else:
        assert list(cache_dir.iterdir()) == []


def test_cache_only_mode_rejects_corrupt_cache_without_network_or_mutation(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
    )
    cache_dir = tmp_path / "cache"
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )
    (cache_dir / "assets" / "js" / "app.js").write_bytes(b"corrupt")
    before = _cache_snapshot(cache_dir)
    fetch = mock.Mock(side_effect=AssertionError("cache-only mode used the network"))
    monkeypatch.setattr(web_cache, "_fetch_web_response", fetch)

    with pytest.raises(ValueError, match="sha256|byte length"):
        web_cache.ensure_web_ui_cached(
            cache_dir=cache_dir,
            source_url="https://viewer.example",
            force=False,
            show_progress=False,
        )

    fetch.assert_not_called()
    assert _cache_snapshot(cache_dir) == before


def test_cache_only_mode_rejects_source_mismatch_without_network_or_mutation(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
    )
    cache_dir = tmp_path / "cache"
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )
    before = _cache_snapshot(cache_dir)
    fetch = mock.Mock(side_effect=AssertionError("cache-only mode used the network"))
    monkeypatch.setattr(web_cache, "_fetch_web_response", fetch)

    with pytest.raises(ValueError, match="Web cache source.*expected"):
        web_cache.ensure_web_ui_cached(
            cache_dir=cache_dir,
            source_url="https://other-viewer.example",
            force=False,
            show_progress=False,
        )

    fetch.assert_not_called()
    assert _cache_snapshot(cache_dir) == before


def test_failed_generation_never_replaces_the_last_complete_cache(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache

    cache_dir = tmp_path / "cache"
    old_assets = _asset_payloads("build-1")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=old_assets,
    )
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )
    before = _cache_snapshot(cache_dir)

    new_assets = _asset_payloads("build-2")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-2",
        assets=new_assets,
        corrupt_path="assets/js/app.js",
    )
    with pytest.raises(ValueError, match="sha256|byte length"):
        web_cache.ensure_web_ui_cached(
            cache_dir=cache_dir,
            source_url="https://viewer.example",
            force=True,
            show_progress=False,
        )

    assert _cache_snapshot(cache_dir) == before
    assert web_cache.verify_web_ui_cache(cache_dir).build_id == "build-1"


@pytest.mark.parametrize(
    ("failure_asset", "expected_attempts"),
    [
        (
            "cellucid-web-assets.json",
            ["cellucid-web-assets.json"],
        ),
        (
            "assets/js/app.js",
            [
                "cellucid-web-assets.json",
                "assets/css/main.css",
                "assets/js/app.js",
            ],
        ),
    ],
)
def test_force_refresh_network_failure_propagates_and_preserves_prior_generation(
    monkeypatch,
    tmp_path: Path,
    failure_asset: str,
    expected_attempts: list[str],
) -> None:
    from cellucid import web_cache

    cache_dir = tmp_path / "cache"
    old_assets = _asset_payloads("build-1")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=old_assets,
    )
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )
    before = _cache_snapshot(cache_dir)

    new_assets = _asset_payloads("build-2")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-2",
        assets=new_assets,
    )
    successful_fetch = web_cache._fetch_web_response
    attempts: list[str] = []

    def fail_selected_request(
        *,
        source_url: str,
        asset_path: str,
        timeout: float,
    ):
        attempts.append(asset_path)
        if asset_path == failure_asset:
            raise ConnectionError("refresh network failure")
        return successful_fetch(
            source_url=source_url,
            asset_path=asset_path,
            timeout=timeout,
        )

    monkeypatch.setattr(web_cache, "_fetch_web_response", fail_selected_request)

    with pytest.raises(ConnectionError, match="refresh network failure"):
        web_cache.ensure_web_ui_cached(
            cache_dir=cache_dir,
            source_url="https://viewer.example",
            force=True,
            show_progress=False,
        )

    assert attempts == expected_attempts
    assert _cache_snapshot(cache_dir) == before
    assert web_cache.verify_web_ui_cache(cache_dir).build_id == "build-1"
    assert list(tmp_path.glob(".cache.staging-*")) == []


def test_cached_corruption_is_rejected_and_never_fetched_on_request(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache
    from cellucid._server_base import CORSMixin

    assets = _asset_payloads("build-1")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
    )
    cache_dir = tmp_path / "cache"
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )

    class Handler(CORSMixin):
        def __init__(self) -> None:
            self.serve_web_ui = True
            self.web_cache_dir = cache_dir
            self.status: int | None = None
            self.headers: dict[str, str] = {}
            self.writer = io.BytesIO()

        def send_response(self, code, message=None) -> None:
            del message
            self.status = int(code)

        def send_header(self, keyword: str, value: object) -> None:
            self.headers[keyword] = str(value)

        def end_headers(self) -> None:
            return None

        def _response_writer(self):
            return self.writer

    with mock.patch(
        "urllib.request.urlopen",
        side_effect=AssertionError("request-time network access"),
    ):
        handler = Handler()
        assert handler.handle_web_asset_get("/assets/js/app.js") is True
        assert handler.status == 200
        assert handler.headers["Content-Type"] == "application/javascript; charset=utf-8"
        assert handler.writer.getvalue() == assets["assets/js/app.js"][0]

        (cache_dir / "assets" / "js" / "app.js").write_bytes(b"changed")
        with pytest.raises(ValueError, match="sha256|byte length"):
            Handler().handle_web_asset_get("/assets/js/app.js")

        assert Handler().handle_web_asset_get("/cellucid-web-assets.json") is False


def test_cache_verification_rejects_symlinked_asset(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
    )
    cache_dir = tmp_path / "cache"
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )

    asset_path = cache_dir / "assets" / "js" / "app.js"
    outside = tmp_path / "outside.js"
    outside.write_bytes(asset_path.read_bytes())
    asset_path.unlink()
    try:
        asset_path.symlink_to(outside)
    except OSError as error:
        pytest.skip(f"Symlinks are unavailable on this platform: {error}")

    with pytest.raises(ValueError, match="symbolic link"):
        web_cache.verify_web_ui_cache(cache_dir)


def test_source_url_must_already_be_canonical(tmp_path: Path) -> None:
    from cellucid.web_cache import ensure_web_ui_cached

    with pytest.raises(ValueError, match="source_url"):
        ensure_web_ui_cached(
            cache_dir=tmp_path / "cache",
            source_url="https://viewer.example/",
            show_progress=False,
        )


def test_cache_establishment_is_the_exact_current_default() -> None:
    import inspect

    from cellucid.jupyter import BaseViewer
    from cellucid.web_cache import ensure_web_ui_cached

    assert inspect.signature(ensure_web_ui_cached).parameters["force"].default is True
    assert (
        inspect.signature(BaseViewer.ensure_web_ui_cached)
        .parameters["force"]
        .default
        is True
    )


def test_server_start_propagates_cache_initialization_failure(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid.server import CellucidServer

    _write_current_prepared_dataset(
        tmp_path,
        dataset_id="cache-failure",
        dataset_name="Cache failure",
    )

    def fail(**kwargs):
        assert kwargs["force"] is True
        raise OSError("cache directory is read-only")

    monkeypatch.setattr("cellucid.web_cache.ensure_web_ui_cached", fail)
    server = CellucidServer(
        tmp_path,
        port=0,
        quiet=True,
        serve_web_ui=True,
        web_source_url="https://viewer.example",
        web_cache_dir=tmp_path / "web-cache",
    )

    with pytest.raises(OSError, match="read-only"):
        server.start_background()
    assert server._server is None
    assert server.is_running() is False


def test_prepared_server_reports_exact_viewer_generation_establishment(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid.server import CellucidServer

    dataset = tmp_path / "dataset"
    _write_current_prepared_dataset(
        dataset,
        dataset_id="generation-report",
        dataset_name="Generation report",
    )
    detail_calls: list[tuple[str, str]] = []
    success_calls: list[str] = []
    cache_calls: list[dict[str, object]] = []

    monkeypatch.setattr(
        "cellucid.server.print_detail",
        lambda label, value: detail_calls.append((label, value)),
    )
    monkeypatch.setattr(
        "cellucid.server.print_success",
        lambda message: success_calls.append(message),
    )
    monkeypatch.setattr(
        "cellucid.web_cache.ensure_web_ui_cached",
        lambda **kwargs: cache_calls.append(kwargs),
    )
    server = CellucidServer(
        dataset,
        host="127.0.0.1",
        port=0,
        quiet=False,
        serve_web_ui=True,
        web_source_url="https://viewer.example",
        web_cache_dir=tmp_path / "web-cache",
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


def test_exported_server_reserves_inventory_and_unlisted_asset_paths(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache
    from cellucid.server import CellucidServer

    assets = _asset_payloads("build-1")
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
    )
    cache_dir = tmp_path / "cache"
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )

    dataset = tmp_path / "dataset"
    (dataset / "assets").mkdir(parents=True)
    (dataset / "assets" / "not-listed.js").write_bytes(b"must not be served")
    (dataset / "cellucid-web-assets.json").write_bytes(b"must not be served")
    _write_current_prepared_dataset(
        dataset,
        dataset_id="reserved-routes",
        dataset_name="Reserved routes",
    )

    server = CellucidServer(
        dataset,
        host="127.0.0.1",
        port=0,
        quiet=True,
        serve_web_ui=True,
        web_source_url="https://viewer.example",
        web_cache_dir=cache_dir,
    )
    server.start_background()
    assert server._server is not None
    host, port = server._server.server_address
    try:
        connection = http.client.HTTPConnection(host, port, timeout=5)
        try:
            connection.request(
                "GET",
                "/?jupyter=true&viewerId=viewer-1&viewerToken=secret-1",
            )
            response = connection.getresponse()
            assert response.status == 200
            assert response.read() == assets["index.html"][0]
        finally:
            connection.close()

        for path in ("/assets/not-listed.js", "/cellucid-web-assets.json"):
            connection = http.client.HTTPConnection(host, port, timeout=5)
            try:
                connection.request("GET", path)
                response = connection.getresponse()
                assert response.status == 404
                response.read()
            finally:
                connection.close()
    finally:
        server.stop()


def test_cache_cleanup_failure_is_propagated(monkeypatch, tmp_path: Path) -> None:
    from cellucid.web_cache import clear_web_cache

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / "index.html").write_text("old", encoding="utf-8")

    monkeypatch.setattr(
        "cellucid.web_cache._remove_web_cache_dir",
        mock.Mock(side_effect=OSError("cannot remove cache")),
    )
    with pytest.raises(OSError, match="cannot remove cache"):
        clear_web_cache(cache_dir=cache_dir)
