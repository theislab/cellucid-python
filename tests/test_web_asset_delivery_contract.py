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
    served_content_types: dict[str, str] | None = None,
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
        if served_content_types is not None:
            # Let a test model a source whose header differs from the header the
            # inventory declares, which is what a real host does.
            content_type = served_content_types.get(asset_path, content_type)
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


def _anndata_server(**kwargs: object):
    import anndata as ad
    import numpy as np
    import pandas as pd

    from cellucid.anndata_server import AnnDataServer

    adata = ad.AnnData(
        X=np.ones((2, 1), dtype=np.float32),
        obs=pd.DataFrame(index=["cell-1", "cell-2"]),
        var=pd.DataFrame(index=["gene-1"]),
        obsm={"X_umap": np.zeros((2, 2), dtype=np.float32)},
    )
    return AnnDataServer(
        adata,
        dataset_name="Web assets",
        dataset_id="web-assets",
        port=0,
        **kwargs,
    )


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


def test_prefetch_accepts_a_source_that_omits_the_charset_parameter(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """A host may format its own Content-Type differently than we serve it.

    GitHub Pages serves ``.xml`` as ``application/xml`` while the inventory
    declares ``application/xml; charset=utf-8`` -- the header this package sends
    when it later serves the asset. Requiring the two to be equal made every
    viewer-startup path fail against the real deployment.
    """
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    assets["browserconfig.xml"] = (
        b"<?xml version='1.0'?><browserconfig/>\n",
        "application/xml; charset=utf-8",
    )
    served = {"browserconfig.xml": "application/xml"}
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
        served_content_types=served,
    )

    web_cache.ensure_web_ui_cached(
        cache_dir=tmp_path / "cache",
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )

    payload, declared = assets["browserconfig.xml"]
    assert (tmp_path / "cache" / "browserconfig.xml").read_bytes() == payload
    # The asset is still served with the inventory's own header, not the
    # source's, so relaxing the fetch check does not change delivery.
    _data, content_type = web_cache.read_cached_web_asset(
        tmp_path / "cache",
        "browserconfig.xml",
    )
    assert content_type == declared


def test_prefetch_still_rejects_a_genuinely_different_media_type(
    monkeypatch,
    tmp_path: Path,
) -> None:
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    assets["browserconfig.xml"] = (
        b"<?xml version='1.0'?><browserconfig/>\n",
        "application/xml; charset=utf-8",
    )
    # A captive portal or error page served with HTTP 200 in place of the asset.
    served = {"browserconfig.xml": "text/html; charset=utf-8"}
    _install_fetch_fixture(
        monkeypatch,
        build_id="build-1",
        assets=assets,
        served_content_types=served,
    )

    with pytest.raises(ValueError, match="response media type"):
        web_cache.ensure_web_ui_cached(
            cache_dir=tmp_path / "cache",
            source_url="https://viewer.example",
            force=True,
            show_progress=False,
        )


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


@pytest.mark.parametrize(
    "source_url",
    [
        "http://viewer.example",
        "http://viewer.example:8080",
        "http://www.cellucid.com",
        "http://192.168.1.10:4173",
        "http://0.0.0.0:4173",
        # Neither of these is loopback; only an exact match may be exempt.
        "http://localhost.viewer.example",
        "http://127.0.0.1.viewer.example",
        # Integer-encoded 127.0.0.1: not a recognised literal, so not exempt.
        "http://2130706433",
    ],
)
def test_plaintext_http_web_source_is_refused_for_non_loopback_hosts(
    source_url: str,
    tmp_path: Path,
) -> None:
    """The unpinned inventory carries the hashes, so http off-machine is rewritable.

    Every asset is SHA-256 verified against ``cellucid-web-assets.json``, but
    that inventory is fetched from the same origin with nothing pinning it. Over
    plaintext HTTP an on-path attacker serves its own assets *and* its own
    hashes, verification passes, and the result is executed as the viewer UI.
    """
    from cellucid.web_cache import ensure_web_ui_cached

    with pytest.raises(ValueError, match="must use https"):
        ensure_web_ui_cached(
            cache_dir=tmp_path / "cache",
            source_url=source_url,
            show_progress=False,
        )


@pytest.mark.parametrize(
    "source_url",
    [
        "http://127.0.0.1:4173",
        "http://127.0.0.1",
        "http://127.1.2.3:9000",
        "http://localhost:4173",
        "http://LOCALHOST:4173",
        "http://[::1]:4173",
        "https://viewer.example",
        "https://www.cellucid.com",
    ],
)
def test_loopback_http_and_any_https_web_source_stay_accepted(source_url: str) -> None:
    """Serving a locally built web app over loopback is a documented path.

    ``docs/user_guide/web_app/b_data_loading/04_server_tutorial.md`` pins
    ``--web-source-url http://127.0.0.1:4173`` for a build from the checked-out
    repository. Loopback traffic never leaves the machine, so there is no
    on-path position for an attacker to occupy.
    """
    from cellucid.web_cache import _require_source_url

    assert _require_source_url(source_url) == source_url


def test_a_cache_recorded_against_a_plaintext_http_source_is_refused(
    monkeypatch,
    tmp_path: Path,
) -> None:
    """The scheme rule also guards the cached source-control file on read."""
    import json as _json

    from cellucid import web_cache

    _install_fetch_fixture(monkeypatch, build_id="build-1", assets=_asset_payloads("build-1"))
    cache_dir = tmp_path / "cache"
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,
        source_url="https://viewer.example",
        force=True,
        show_progress=False,
    )

    (cache_dir / ".cellucid-web-source.json").write_text(
        _json.dumps({"version": 1, "source_url": "http://viewer.example"}, separators=(",", ":")),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="must use https"):
        web_cache.verify_web_ui_cache(cache_dir)


def test_cache_establishment_is_the_exact_current_default() -> None:
    import inspect

    from cellucid.jupyter import AnnDataViewer
    from cellucid.web_cache import ensure_web_ui_cached

    assert inspect.signature(ensure_web_ui_cached).parameters["force"].default is True
    assert (
        inspect.signature(AnnDataViewer.ensure_web_ui_cached)
        .parameters["force"]
        .default
        is True
    )


def test_server_start_propagates_cache_initialization_failure(
    monkeypatch,
    tmp_path: Path,
) -> None:
    def fail(**kwargs):
        assert kwargs["force"] is True
        raise OSError("cache directory is read-only")

    monkeypatch.setattr("cellucid.web_cache.ensure_web_ui_cached", fail)
    server = _anndata_server(
        quiet=True,
        serve_web_ui=True,
        web_source_url="https://viewer.example",
        web_cache_dir=tmp_path / "web-cache",
    )

    with pytest.raises(OSError, match="read-only"):
        server.start_background()
    assert server._server is None
    assert server.is_running() is False


def test_the_server_reports_exact_viewer_generation_establishment(
    monkeypatch,
    tmp_path: Path,
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
    server = _anndata_server(
        host="127.0.0.1",
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


def test_the_server_reserves_inventory_and_unlisted_asset_paths(
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

    server = _anndata_server(
        host="127.0.0.1",
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
                "/?viewer=true&viewerId=viewer-1&viewerToken=secret-1",
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


def test_an_asset_is_verified_once_per_version_not_on_every_request(monkeypatch, tmp_path: Path) -> None:
    """Hashing every asset on every request cost 30 ms each: seconds before a first frame."""
    from cellucid import web_cache

    assets = _asset_payloads("build-1")
    _install_fetch_fixture(monkeypatch, build_id="build-1", assets=assets)
    cache_dir = tmp_path / "cache"
    web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir, source_url="https://viewer.example", force=True, show_progress=False
    )
    hashed: list[str] = []
    loaded: list[Path] = []
    verify, load = web_cache._verify_asset_bytes, web_cache._load_cached_inventory
    monkeypatch.setattr(web_cache, "_verify_asset_bytes", lambda asset, data: (hashed.append(asset.path), verify(asset, data)))
    monkeypatch.setattr(web_cache, "_load_cached_inventory", lambda path, **kw: (loaded.append(path), load(path, **kw))[1])
    for _ in range(5):
        assert web_cache.read_cached_web_asset(cache_dir, "assets/js/app.js")[0] == assets["assets/js/app.js"][0]
    assert hashed == ["assets/js/app.js"]
    assert len(loaded) == 1
    # A replaced file is a new version: checked again, and refused when wrong.
    target = cache_dir / "assets" / "js" / "app.js"
    target.unlink()
    target.write_bytes(b"x" * len(assets["assets/js/app.js"][0]))
    with pytest.raises(ValueError, match="sha256"):
        web_cache.read_cached_web_asset(cache_dir, "assets/js/app.js")
