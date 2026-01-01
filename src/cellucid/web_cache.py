"""
Cellucid web UI cache utilities.

Cellucid notebooks embed the viewer UI from the local Python data server. When
the web UI is not packaged locally, the server runs in "hosted-asset proxy"
mode: it downloads `index.html` + `/assets/*` from `CELLUCID_WEB_URL` and caches
those files on disk.

This module exposes:
- `get_web_cache_dir()`: where cached assets live
- `clear_web_cache()`: delete cached assets (forces a re-download)
- `ensure_web_ui_cached()`: best-effort prefetch of the full UI asset graph
"""

from __future__ import annotations

import json
import logging
import posixpath
import re
import time
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ._server_base import (
    CELLUCID_WEB_URL,
    _extract_web_build_id,
    _is_proxyable_web_path,
    _normalize_web_path,
    _purge_web_proxy_cache_dir,
    _web_proxy_cache_dir,
)

logger = logging.getLogger("cellucid.web_cache")

_PREFETCH_META_NAME = ".cellucid-web-prefetch.json"
_SCAN_TEXT_MAX_BYTES = 8 * 1024 * 1024  # 8MB cap for heuristic string scanning


@dataclass(frozen=True)
class WebCachePrefetchSummary:
    cache_dir: Path
    source_url: str
    build_id: str | None
    downloaded_files: int
    downloaded_bytes: int
    skipped_files: int
    errors: list[str]


def get_web_cache_dir() -> Path:
    """Return the on-disk cache directory used for hosted-asset proxy mode."""
    return _web_proxy_cache_dir()


def clear_web_cache(*, cache_dir: Path | None = None) -> Path:
    """
    Delete cached web UI assets.

    Returns the cache directory path (even if it didn't exist).
    """
    cache_dir = (cache_dir or _web_proxy_cache_dir()).expanduser().resolve()
    _purge_web_proxy_cache_dir(cache_dir)
    return cache_dir


def _prefetch_meta_path(cache_dir: Path) -> Path:
    return cache_dir / _PREFETCH_META_NAME


def _read_prefetch_meta(cache_dir: Path) -> dict[str, Any] | None:
    path = _prefetch_meta_path(cache_dir)
    try:
        if not path.is_file():
            return None
        data = json.loads(path.read_text("utf-8"))
        return data if isinstance(data, dict) else None
    except Exception:
        return None


def _write_prefetch_meta(cache_dir: Path, meta: dict[str, Any]) -> None:
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        _prefetch_meta_path(cache_dir).write_text(json.dumps(meta, indent=2, sort_keys=True), "utf-8")
    except Exception:
        # Cache is best-effort.
        pass


def _extract_candidate_paths(text: str) -> set[str]:
    """
    Deprecated wrapper for asset discovery without file context.

    Prefer `_extract_candidate_paths_for_file(text, base_url_path=...)` so
    relative imports (e.g. `../data/foo.js`) can be resolved correctly.
    """
    return _extract_candidate_paths_for_file(text, base_url_path="/")


def _resolve_relative_web_path(raw: str, *, base_url_path: str) -> str | None:
    raw = (raw or "").strip()
    if not raw:
        return None

    # Drop query/hash for cache purposes.
    raw = raw.split("#", 1)[0].split("?", 1)[0].strip()
    if not raw:
        return None

    # Ignore obvious non-path schemes.
    lowered = raw.lower()
    if lowered.startswith(("http://", "https://", "data:", "blob:", "ws:", "wss:", "mailto:", "tel:")):
        return None

    if raw.startswith("//"):
        return None

    # Root-relative.
    if raw.startswith("/"):
        return raw

    # Common build output: "assets/..." is meant to be root-relative.
    if raw.startswith("assets/"):
        return "/" + raw.lstrip("/")

    # Resolve relative to the file's directory.
    base = base_url_path or "/"
    base_dir = base if base.endswith("/") else (posixpath.dirname(base) + "/")
    resolved = posixpath.normpath(posixpath.join(base_dir, raw))
    if not resolved.startswith("/"):
        resolved = "/" + resolved
    return resolved


def _strip_js_comments(text: str) -> str:
    """
    Best-effort JavaScript comment stripper used to avoid prefetching paths from
    JSDoc/type comments (e.g. `import('./foo.js')`).
    """
    out: list[str] = []
    i = 0
    n = len(text)
    state = "code"  # code | sq | dq | bt

    while i < n:
        ch = text[i]
        nxt = text[i + 1] if i + 1 < n else ""

        if state == "code":
            if ch == "'" and nxt:
                state = "sq"
                out.append(ch)
                i += 1
                continue
            if ch == '"' and nxt:
                state = "dq"
                out.append(ch)
                i += 1
                continue
            if ch == "`":
                state = "bt"
                out.append(ch)
                i += 1
                continue
            if ch == "/" and nxt == "/":
                # Line comment: skip to newline (but keep the newline).
                i += 2
                while i < n and text[i] not in "\r\n":
                    i += 1
                continue
            if ch == "/" and nxt == "*":
                # Block comment: skip to closing.
                i += 2
                while i + 1 < n and not (text[i] == "*" and text[i + 1] == "/"):
                    i += 1
                i += 2 if i + 1 < n else 0
                continue
            out.append(ch)
            i += 1
            continue

        if state == "sq":
            out.append(ch)
            if ch == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            if ch == "'":
                state = "code"
            i += 1
            continue

        if state == "dq":
            out.append(ch)
            if ch == "\\" and i + 1 < n:
                out.append(text[i + 1])
                i += 2
                continue
            if ch == '"':
                state = "code"
            i += 1
            continue

        # state == "bt"
        out.append(ch)
        if ch == "\\" and i + 1 < n:
            out.append(text[i + 1])
            i += 2
            continue
        if ch == "`":
            state = "code"
        i += 1

    return "".join(out)


def _extract_candidate_paths_for_file(text: str, *, base_url_path: str) -> set[str]:
    out: set[str] = set()
    root_files = {
        "/index.html",
        "/favicon.ico",
        "/site.webmanifest",
        "/browserconfig.xml",
    }
    allowed_exts = (
        ".html",
        ".js",
        ".css",
        ".json",
        ".svg",
        ".png",
        ".jpg",
        ".jpeg",
        ".ico",
        ".xml",
        ".webmanifest",
        ".woff2",
        ".woff",
        ".ttf",
        ".otf",
        ".wasm",
        ".bin",
        ".gz",
        ".txt",
        ".map",
    )

    def add(raw: str) -> None:
        resolved = _resolve_relative_web_path(raw, base_url_path=base_url_path)
        if resolved:
            if resolved not in root_files:
                low = resolved.lower()
                if not low.endswith(allowed_exts):
                    return
            out.add(resolved)

    ext = posixpath.splitext(base_url_path or "")[1].lower()

    if ext in (".html", ".htm"):
        for m in re.finditer(r"""(?:(?:src|href)=["'])([^"'?#]+)""", text, flags=re.IGNORECASE):
            add(m.group(1))

    if ext == ".css":
        for m in re.finditer(r"""url\(\s*["']?([^"')?#]+)""", text):
            add(m.group(1))
        for m in re.finditer(r"""@import\s+["']([^"']+)""", text, flags=re.IGNORECASE):
            add(m.group(1))

    if ext in (".js", ".mjs"):
        code = _strip_js_comments(text)
        for m in re.finditer(r"""(?m)\bfrom\s+["']([^"']+)["']""", code):
            add(m.group(1))
        for m in re.finditer(r"""(?m)\bimport\s+["']([^"']+)["']""", code):
            add(m.group(1))
        for m in re.finditer(r"""\bimport\s*\(\s*["']([^"']+)["']\s*\)""", code):
            add(m.group(1))
        for m in re.finditer(r"""new\s+URL\s*\(\s*["']([^"']+)["']\s*,\s*import\.meta\.url\s*\)""", code):
            add(m.group(1))

    # Conservative fallback: literal /assets/... and assets/... occurrences.
    for m in re.finditer(r"""(/assets/[^"'\s)]+)""", text):
        add(m.group(1))
    for m in re.finditer(r"""\bassets/[^"'\s)]+""", text):
        add(m.group(0))

    # Root files that are commonly referenced but not under /assets/.
    out.update(root_files)

    return out


def _cache_path_for_web_path(cache_dir: Path, url_path: str) -> Path:
    rel = url_path.lstrip("/")
    candidate = (cache_dir / rel).resolve()
    candidate.relative_to(cache_dir.resolve())
    return candidate


def _fetch_to_cache(
    *,
    source_url: str,
    cache_dir: Path,
    url_path: str,
    timeout: float = 15.0,
    user_agent: str = "cellucid-python (web-cache-prefetch)",
) -> tuple[int, str | None]:
    import urllib.request

    remote_url = f"{source_url}{url_path}"
    req = urllib.request.Request(remote_url, headers={"User-Agent": user_agent})

    cache_path = _cache_path_for_web_path(cache_dir, url_path)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    with urllib.request.urlopen(req, timeout=timeout) as resp:
        status = getattr(resp, "status", 200)
        try:
            status_int = 200 if status is None else int(status)
        except Exception:
            status_int = 200
        if status_int != 200:
            raise RuntimeError(f"HTTP {status} for {remote_url}")

        content_type = resp.headers.get("Content-Type")
        total = 0
        with cache_path.open("wb") as f:
            while True:
                chunk = resp.read(64 * 1024)
                if not chunk:
                    break
                f.write(chunk)
                total += len(chunk)

    return total, content_type


def ensure_web_ui_cached(
    *,
    cache_dir: Path | None = None,
    source_url: str = CELLUCID_WEB_URL,
    force: bool = False,
    show_progress: bool = True,
) -> WebCachePrefetchSummary:
    """
    Best-effort: download and cache *all* viewer UI assets for offline-safe use.

    This function is intentionally conservative and dev-phase:
    - It discovers assets by fetching `index.html`, then scanning HTML/JS/CSS for
      referenced files (including module imports and `url(...)` references) and
      recursively downloading them.
    - It does not attempt directory listings.

    Parameters
    ----------
    cache_dir
        Cache directory (defaults to `CELLUCID_WEB_PROXY_CACHE_DIR` or system temp).
    source_url
        Hosted web origin to prefetch from (defaults to `CELLUCID_WEB_URL`).
    force
        If True, re-download even if a prefetch marker exists for this build id.
    show_progress
        If True, display a `tqdm` progress bar (works in terminals and notebooks).
    """
    cache_dir = (cache_dir or _web_proxy_cache_dir()).expanduser().resolve()
    source_url = (source_url or "").rstrip("/")
    errors: list[str] = []

    # Read cached build id (if present).
    cached_index_path = cache_dir / "index.html"
    cached_build_id: str | None = None
    if cached_index_path.is_file():
        try:
            cached_build_id = _extract_web_build_id(cached_index_path.read_bytes())
        except Exception:
            cached_build_id = None

    # Fetch the latest index.html to detect build changes + seed the graph.
    remote_build_id: str | None = None
    try:
        bytes_written, _ct = _fetch_to_cache(
            source_url=source_url,
            cache_dir=cache_dir,
            url_path="/index.html",
            timeout=15.0,
        )
        del bytes_written, _ct
        try:
            remote_build_id = _extract_web_build_id((cache_dir / "index.html").read_bytes())
        except Exception:
            remote_build_id = None
    except Exception as e:
        # Offline path: if we have a cached index.html, proceed with it.
        if not cached_index_path.is_file():
            raise
        errors.append(f"Failed to refresh hosted index.html; using cached copy: {e}")
        remote_build_id = cached_build_id

    # Invalidate cache if the hosted build changed (same rule as the runtime proxy).
    if remote_build_id and cached_build_id and remote_build_id != cached_build_id:
        _purge_web_proxy_cache_dir(cache_dir)
        # Re-fetch index.html after purge.
        try:
            _fetch_to_cache(source_url=source_url, cache_dir=cache_dir, url_path="/index.html", timeout=15.0)
        except Exception as e:
            errors.append(f"Failed to re-download index.html after cache purge: {e}")

    build_id = remote_build_id or cached_build_id

    # Skip if we already prefetched this build.
    meta = _read_prefetch_meta(cache_dir)
    if not force and isinstance(meta, dict):
        meta_build = meta.get("build_id")
        meta_source = meta.get("source_url")
        meta_errors = meta.get("errors")
        meta_ok = meta_errors in (None, [], "")
        if meta_build and meta_build == build_id and meta_source == source_url and meta_ok:
            return WebCachePrefetchSummary(
                cache_dir=cache_dir,
                source_url=source_url,
                build_id=build_id,
                downloaded_files=0,
                downloaded_bytes=0,
                skipped_files=0,
                errors=errors,
            )

    # Seed discovery from cached index.html (now present).
    try:
        index_text = (cache_dir / "index.html").read_text("utf-8", errors="ignore")
    except Exception:
        index_text = ""

    seen: set[str] = set()
    q: deque[str] = deque()

    for p in _extract_candidate_paths_for_file(index_text, base_url_path="/index.html"):
        norm = _normalize_web_path(p)
        if not norm:
            continue
        q.append(norm)

    # Always include index + root files even if not referenced explicitly.
    for p in ("/index.html", "/favicon.ico", "/site.webmanifest", "/browserconfig.xml"):
        q.append(p)

    downloaded_files = 0
    downloaded_bytes = 0
    skipped_files = 0

    bar = None
    if show_progress:
        try:
            from tqdm.auto import tqdm  # type: ignore

            bar = tqdm(
                total=None,
                desc="Cellucid web UI cache",
                unit="file",
            )
        except Exception:
            bar = None

    while q:
        url_path = q.popleft()
        if url_path in seen:
            continue
        seen.add(url_path)

        norm = _normalize_web_path(url_path)
        if not norm:
            continue
        if not _is_proxyable_web_path(norm):
            continue

        cache_path: Path | None = None
        try:
            cache_path = _cache_path_for_web_path(cache_dir, norm)
        except Exception:
            continue

        if cache_path.is_file() and not force:
            skipped_files += 1
            # Still scan cached text files to discover additional lazy chunks
            # that may not have been cached previously.
            try:
                is_textish = norm.endswith((".html", ".js", ".css", ".json", ".svg", ".xml", ".txt"))
                if is_textish:
                    data = cache_path.read_bytes()
                    if len(data) <= _SCAN_TEXT_MAX_BYTES:
                        text = data.decode("utf-8", errors="ignore")
                        for p in _extract_candidate_paths_for_file(text, base_url_path=norm):
                            n2 = _normalize_web_path(p)
                            if n2 and n2 not in seen:
                                q.append(n2)
            except Exception:
                pass
            if bar is not None:
                bar.update(1)
            continue

        try:
            bytes_written, content_type = _fetch_to_cache(
                source_url=source_url,
                cache_dir=cache_dir,
                url_path=norm,
                timeout=15.0,
            )
            downloaded_files += 1
            downloaded_bytes += int(bytes_written)
        except Exception as e:
            errors.append(f"Failed to fetch {norm}: {e}")
            if bar is not None:
                bar.update(1)
            continue

        # Best-effort: recursively discover more assets by scanning text files.
        is_textish = False
        if content_type:
            ct = content_type.lower()
            is_textish = any(
                x in ct
                for x in (
                    "text/",
                    "javascript",
                    "json",
                    "xml",
                    "svg",
                    "css",
                )
            )
        if not is_textish:
            # Fall back to extension checks.
            is_textish = norm.endswith((".html", ".js", ".css", ".json", ".svg", ".xml", ".txt"))

        if is_textish and cache_path is not None:
            try:
                data = cache_path.read_bytes()
                if len(data) <= _SCAN_TEXT_MAX_BYTES:
                    text = data.decode("utf-8", errors="ignore")
                    for p in _extract_candidate_paths_for_file(text, base_url_path=norm):
                        n2 = _normalize_web_path(p)
                        if n2 and n2 not in seen:
                            q.append(n2)
            except Exception:
                pass

        if bar is not None:
            bar.update(1)

    if bar is not None:
        try:
            bar.close()
        except Exception:
            pass

    _write_prefetch_meta(
        cache_dir,
        {
            "source_url": source_url,
            "build_id": build_id,
            "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "downloaded_files": downloaded_files,
            "downloaded_bytes": downloaded_bytes,
            "skipped_files": skipped_files,
            "errors": errors,
        },
    )

    return WebCachePrefetchSummary(
        cache_dir=cache_dir,
        source_url=source_url,
        build_id=build_id,
        downloaded_files=downloaded_files,
        downloaded_bytes=downloaded_bytes,
        skipped_files=skipped_files,
        errors=errors,
    )
