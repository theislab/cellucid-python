from __future__ import annotations

import json
from pathlib import Path


def test_web_cache_prefetch_retries_when_meta_has_errors(monkeypatch, tmp_path: Path) -> None:
    """
    If a previous prefetch recorded errors, we should not treat the cache as complete
    and should attempt to fetch missing assets again.
    """
    from cellucid import web_cache

    source_url = "https://example.com"
    index_html = (
        b"<!doctype html><html><head>"
        b'<meta name="cellucid-web-build-id" content="build-1" />'
        b"</head><body>"
        b'<script type="module" src="assets/js/app/main.js"></script>'
        b"</body></html>"
    )

    # Seed a cached index.html so offline paths still work.
    (tmp_path / "index.html").write_bytes(index_html)

    # Simulate a previous run that had errors (incomplete prefetch).
    (tmp_path / ".cellucid-web-prefetch.json").write_text(
        json.dumps({"source_url": source_url, "build_id": "build-1", "errors": ["timeout"]}),
        "utf-8",
    )

    calls: list[str] = []

    def fake_fetch_to_cache(*, source_url: str, cache_dir: Path, url_path: str, timeout: float = 15.0, user_agent: str = ""):
        calls.append(url_path)
        out_path = cache_dir / url_path.lstrip("/")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if url_path == "/index.html":
            out_path.write_bytes(index_html)
            return len(index_html), "text/html; charset=utf-8"
        if url_path.endswith(".js"):
            payload = b"console.log('ok')\n"
            out_path.write_bytes(payload)
            return len(payload), "application/javascript"
        out_path.write_bytes(b"")
        return 0, "application/octet-stream"

    monkeypatch.setattr(web_cache, "_fetch_to_cache", fake_fetch_to_cache)

    web_cache.ensure_web_ui_cached(
        cache_dir=tmp_path,
        source_url=source_url,
        force=False,
        show_progress=False,
    )

    # The function always refreshes /index.html; the important part is that it
    # does NOT return early and fetches at least one additional asset.
    assert "/index.html" in calls
    assert any(p != "/index.html" for p in calls)


def test_web_cache_prefetch_skips_when_meta_ok(monkeypatch, tmp_path: Path) -> None:
    """
    If the prefetch meta matches the current build and has no errors, we should
    skip the expensive graph walk and only refresh index.html for build checks.
    """
    from cellucid import web_cache

    source_url = "https://example.com"
    index_html = (
        b"<!doctype html><html><head>"
        b'<meta name="cellucid-web-build-id" content="build-1" />'
        b"</head><body>"
        b'<script type="module" src="assets/js/app/main.js"></script>'
        b"</body></html>"
    )
    (tmp_path / "index.html").write_bytes(index_html)

    (tmp_path / ".cellucid-web-prefetch.json").write_text(
        json.dumps({"source_url": source_url, "build_id": "build-1", "errors": []}),
        "utf-8",
    )

    calls: list[str] = []

    def fake_fetch_to_cache(*, source_url: str, cache_dir: Path, url_path: str, timeout: float = 15.0, user_agent: str = ""):
        calls.append(url_path)
        out_path = cache_dir / url_path.lstrip("/")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_bytes(index_html)
        return len(index_html), "text/html; charset=utf-8"

    monkeypatch.setattr(web_cache, "_fetch_to_cache", fake_fetch_to_cache)

    web_cache.ensure_web_ui_cached(
        cache_dir=tmp_path,
        source_url=source_url,
        force=False,
        show_progress=False,
    )

    # Only index refresh should occur (no full traversal when meta is clean).
    assert calls == ["/index.html"]

