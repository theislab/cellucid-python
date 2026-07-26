"""Exact, generation-atomic delivery of the Cellucid web application."""

from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
import threading
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import quote, urlparse

from ._server_base import (
    CELLUCID_WEB_URL,
    WEB_ASSET_INVENTORY_FILENAME,
    _extract_web_build_id,
    _remove_web_cache_dir,
    _require_web_cache_directory_or_missing,
    _web_cache_dir,
)

_WEB_CACHE_SOURCE_FILENAME = ".cellucid-web-source.json"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_MIME_RE = re.compile(
    r"^[a-z0-9][a-z0-9!#$&^_.+-]*/[a-z0-9][a-z0-9!#$&^_.+-]*"
    r"(?:; [a-z0-9][a-z0-9!#$&^_.+-]*=[a-z0-9][a-z0-9!#$&^_.+-]*)*$"
)
_WINDOWS_RESERVED_NAMES = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{number}" for number in range(1, 10)),
    *(f"LPT{number}" for number in range(1, 10)),
}
_CACHE_LOCK = threading.RLock()


@dataclass(frozen=True)
class WebAsset:
    """One immutable object declared by the web build inventory."""

    path: str
    sha256: str
    bytes: int
    content_type: str


@dataclass(frozen=True)
class WebAssetInventory:
    """The exact published web build inventory."""

    version: int
    build_id: str
    assets: tuple[WebAsset, ...]


@dataclass(frozen=True)
class FetchedWebResponse:
    """A fully read HTTP response with required delivery metadata."""

    data: bytes
    content_type: str
    content_length: int


@dataclass(frozen=True)
class WebCachePrefetchSummary:
    """Result of establishing one complete local web build."""

    cache_dir: Path
    source_url: str
    build_id: str
    downloaded_files: int
    downloaded_bytes: int
    skipped_files: int


def _require_source_url(value: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
        or any(character.isspace() for character in value)
    ):
        raise ValueError("source_url must be a non-empty URL without surrounding whitespace")
    if value.endswith("/"):
        raise ValueError("source_url must not end with '/'")
    parsed = urlparse(value)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise ValueError("source_url must be an absolute HTTP or HTTPS URL")
    if parsed.username is not None or parsed.password is not None:
        raise ValueError("source_url must not contain credentials")
    if parsed.query or parsed.fragment:
        raise ValueError("source_url must not contain a query or fragment")
    try:
        if parsed.hostname is None:
            raise ValueError("source_url must contain a hostname")
        _ = parsed.port
    except ValueError as error:
        raise ValueError("source_url contains an invalid host or port") from error
    return value


def _require_cache_dir(value: str | Path | None) -> Path:
    cache_dir = Path(value) if value is not None else _web_cache_dir()
    cache_dir = Path(os.path.abspath(cache_dir.expanduser()))
    if cache_dir == Path(cache_dir.anchor):
        raise ValueError("The filesystem root cannot be used as the web cache directory")
    _require_web_cache_directory_or_missing(cache_dir)
    return cache_dir


def _require_content_type(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or not _MIME_RE.fullmatch(value):
        raise ValueError(f"{label} must be an exact normalized MIME type")
    return value


def _require_asset_path(value: Any) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError("Web asset path must be a non-empty string")
    if value.startswith("/") or "\\" in value or "?" in value or "#" in value:
        raise ValueError(f"Web asset path is not root-relative POSIX: {value!r}")
    if any(character.isspace() for character in value):
        raise ValueError(f"Web asset path contains whitespace: {value!r}")
    if any(ord(character) < 32 or ord(character) == 127 for character in value):
        raise ValueError(f"Web asset path contains a control character: {value!r}")
    if any(character in '<>:"|*' for character in value):
        raise ValueError(f"Web asset path is not portable across platforms: {value!r}")
    parts = value.split("/")
    if any(part in {"", ".", ".."} for part in parts):
        raise ValueError(f"Web asset path contains an empty or dot segment: {value!r}")
    for part in parts:
        if part.endswith(".") or part.split(".", 1)[0].upper() in _WINDOWS_RESERVED_NAMES:
            raise ValueError(
                f"Web asset path is not portable across platforms: {value!r}"
            )
    if value == WEB_ASSET_INVENTORY_FILENAME:
        raise ValueError("The web asset inventory cannot declare itself as an asset")
    return value


def parse_web_asset_inventory(payload: bytes) -> WebAssetInventory:
    """Parse the one accepted web asset inventory schema."""
    if not isinstance(payload, bytes):
        raise TypeError("Web asset inventory payload must be bytes")
    try:
        raw = json.loads(payload.decode("utf-8"))
    except UnicodeDecodeError as error:
        raise ValueError("Web asset inventory must be valid UTF-8") from error
    except json.JSONDecodeError as error:
        raise ValueError("Web asset inventory must be valid JSON") from error

    if not isinstance(raw, dict) or set(raw) != {"version", "build_id", "assets"}:
        raise ValueError(
            "Web asset inventory must contain exactly version, build_id, and assets"
        )
    if type(raw["version"]) is not int or raw["version"] != 1:
        raise ValueError("Web asset inventory version must be integer 1")
    build_id = raw["build_id"]
    if (
        not isinstance(build_id, str)
        or not build_id
        or build_id != build_id.strip()
        or len(build_id) > 200
    ):
        raise ValueError("Web asset inventory build_id must be a non-empty exact string")
    raw_assets = raw["assets"]
    if not isinstance(raw_assets, list) or not raw_assets:
        raise ValueError("Web asset inventory assets must be a non-empty array")

    assets: list[WebAsset] = []
    for index, raw_asset in enumerate(raw_assets):
        if not isinstance(raw_asset, dict) or set(raw_asset) != {
            "path",
            "sha256",
            "bytes",
            "content_type",
        }:
            raise ValueError(
                f"Web asset inventory entry {index} must contain exactly "
                "path, sha256, bytes, and content_type"
            )
        path = _require_asset_path(raw_asset["path"])
        sha256 = raw_asset["sha256"]
        if not isinstance(sha256, str) or not _SHA256_RE.fullmatch(sha256):
            raise ValueError(f"Web asset {path!r} has an invalid sha256")
        byte_count = raw_asset["bytes"]
        if type(byte_count) is not int or byte_count < 0:
            raise ValueError(f"Web asset {path!r} has an invalid byte length")
        content_type = _require_content_type(
            raw_asset["content_type"],
            label=f"Web asset {path!r} content_type",
        )
        assets.append(
            WebAsset(
                path=path,
                sha256=sha256,
                bytes=byte_count,
                content_type=content_type,
            )
        )

    paths = [asset.path for asset in assets]
    if paths != sorted(paths):
        raise ValueError("Web asset inventory paths must be sorted")
    if len(paths) != len(set(paths)):
        raise ValueError("Web asset inventory paths must be unique")
    if "index.html" not in set(paths):
        raise ValueError("Web asset inventory must contain index.html")

    return WebAssetInventory(version=1, build_id=build_id, assets=tuple(assets))


def _source_control_bytes(source_url: str) -> bytes:
    return json.dumps(
        {"version": 1, "source_url": source_url},
        separators=(",", ":"),
    ).encode("utf-8")


def _read_cache_source(cache_dir: Path) -> str:
    source_path = cache_dir / _WEB_CACHE_SOURCE_FILENAME
    try:
        raw = json.loads(source_path.read_text(encoding="utf-8"))
    except FileNotFoundError as error:
        raise ValueError(f"Web cache source control is missing: {source_path}") from error
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"Web cache source control is invalid: {source_path}") from error
    if not isinstance(raw, dict) or set(raw) != {"version", "source_url"}:
        raise ValueError("Web cache source control has a noncanonical schema")
    if type(raw["version"]) is not int or raw["version"] != 1:
        raise ValueError("Web cache source control version must be integer 1")
    return _require_source_url(raw["source_url"])


def _verify_asset_bytes(asset: WebAsset, data: bytes) -> None:
    if len(data) != asset.bytes:
        raise ValueError(
            f"Web asset {asset.path!r} byte length is {len(data)}, expected {asset.bytes}"
        )
    digest = hashlib.sha256(data).hexdigest()
    if digest != asset.sha256:
        raise ValueError(
            f"Web asset {asset.path!r} sha256 is {digest}, expected {asset.sha256}"
        )


def _load_cached_inventory(
    cache_dir: Path,
    *,
    expected_source_url: str | None = None,
) -> WebAssetInventory:
    if not cache_dir.is_dir():
        raise ValueError(f"Web cache directory does not exist: {cache_dir}")
    source_url = _read_cache_source(cache_dir)
    if expected_source_url is not None:
        expected = _require_source_url(expected_source_url)
        if source_url != expected:
            raise ValueError(
                f"Web cache source is {source_url!r}, expected {expected!r}"
            )
    inventory_path = cache_dir / WEB_ASSET_INVENTORY_FILENAME
    try:
        inventory = parse_web_asset_inventory(inventory_path.read_bytes())
    except FileNotFoundError as error:
        raise ValueError(f"Web asset inventory is missing: {inventory_path}") from error
    index_path = cache_dir / "index.html"
    try:
        index_build_id = _extract_web_build_id(index_path.read_bytes())
    except FileNotFoundError as error:
        raise ValueError(f"Web cache index is missing: {index_path}") from error
    if index_build_id != inventory.build_id:
        raise ValueError(
            f"Web cache index build id is {index_build_id!r}, "
            f"inventory declares {inventory.build_id!r}"
        )
    return inventory


def verify_web_ui_cache(
    cache_dir: str | Path,
    *,
    expected_source_url: str | None = None,
) -> WebAssetInventory:
    """Verify every control file and every declared object in a cache generation."""
    resolved = _require_cache_dir(cache_dir)
    with _CACHE_LOCK:
        inventory = _load_cached_inventory(
            resolved,
            expected_source_url=expected_source_url,
        )
        for path in resolved.rglob("*"):
            if path.is_symlink():
                raise ValueError(f"Web cache contains a symbolic link: {path}")
        expected_files = {
            WEB_ASSET_INVENTORY_FILENAME,
            _WEB_CACHE_SOURCE_FILENAME,
            *(asset.path for asset in inventory.assets),
        }
        actual_files = {
            path.relative_to(resolved).as_posix()
            for path in resolved.rglob("*")
            if path.is_file()
        }
        if actual_files != expected_files:
            missing = sorted(expected_files - actual_files)
            unexpected = sorted(actual_files - expected_files)
            raise ValueError(
                f"Web cache file set differs from its inventory; "
                f"missing={missing}, unexpected={unexpected}"
            )
        for asset in inventory.assets:
            _verify_asset_bytes(asset, (resolved / asset.path).read_bytes())
        return inventory


def read_cached_web_asset(
    cache_dir: str | Path,
    asset_path: str,
) -> tuple[bytes, str] | None:
    """Read and verify one inventory-declared object without network access."""
    resolved = _require_cache_dir(cache_dir)
    normalized_path = _require_asset_path(asset_path)
    with _CACHE_LOCK:
        inventory = _load_cached_inventory(resolved)
        asset_by_path = {asset.path: asset for asset in inventory.assets}
        asset = asset_by_path.get(normalized_path)
        if asset is None:
            return None
        data = (resolved / normalized_path).read_bytes()
        _verify_asset_bytes(asset, data)
        return data, asset.content_type


def _fetch_web_response(
    *,
    source_url: str,
    asset_path: str,
    timeout: float,
) -> FetchedWebResponse:
    import urllib.request

    source = _require_source_url(source_url)
    path = (
        WEB_ASSET_INVENTORY_FILENAME
        if asset_path == WEB_ASSET_INVENTORY_FILENAME
        else _require_asset_path(asset_path)
    )
    if not isinstance(timeout, int | float) or isinstance(timeout, bool) or timeout <= 0:
        raise ValueError("timeout must be a positive number")
    remote_url = f"{source}/{quote(path, safe='/-._~')}"
    request = urllib.request.Request(
        remote_url,
        headers={"User-Agent": "cellucid-python/web-cache"},
    )
    with urllib.request.urlopen(request, timeout=float(timeout)) as response:
        status = response.status
        if type(status) is not int or status != 200:
            raise RuntimeError(f"Expected HTTP 200 for {remote_url}, received {status!r}")
        content_type = _require_content_type(
            response.headers.get("Content-Type"),
            label=f"HTTP Content-Type for {remote_url}",
        )
        raw_content_length = response.headers.get("Content-Length")
        if not isinstance(raw_content_length, str) or not raw_content_length.isdecimal():
            raise ValueError(f"HTTP Content-Length is missing or invalid for {remote_url}")
        content_length = int(raw_content_length)
        data = response.read()
    if len(data) != content_length:
        raise ValueError(
            f"HTTP body length for {remote_url} is {len(data)}, expected {content_length}"
        )
    return FetchedWebResponse(
        data=data,
        content_type=content_type,
        content_length=content_length,
    )


def _require_fetched_response(
    response: FetchedWebResponse,
    *,
    label: str,
) -> FetchedWebResponse:
    if type(response) is not FetchedWebResponse:
        raise TypeError(f"{label} must be returned as FetchedWebResponse")
    if not isinstance(response.data, bytes):
        raise TypeError(f"{label} body must be bytes")
    _require_content_type(response.content_type, label=f"{label} Content-Type")
    if type(response.content_length) is not int or response.content_length < 0:
        raise ValueError(f"{label} Content-Length must be a non-negative integer")
    if len(response.data) != response.content_length:
        raise ValueError(
            f"{label} body length is {len(response.data)}, "
            f"Content-Length is {response.content_length}"
        )
    return response


def _verify_fetched_asset(asset: WebAsset, response: FetchedWebResponse) -> None:
    _require_fetched_response(response, label=f"Web asset {asset.path!r}")
    if response.content_length != asset.bytes:
        raise ValueError(
            f"Web asset {asset.path!r} response byte length is "
            f"{response.content_length}, expected {asset.bytes}"
        )
    if response.content_type != asset.content_type:
        raise ValueError(
            f"Web asset {asset.path!r} response content type is "
            f"{response.content_type!r}, expected {asset.content_type!r}"
        )
    _verify_asset_bytes(asset, response.data)


def _establish_staged_generation(
    *,
    cache_dir: Path,
    source_url: str,
    inventory: WebAssetInventory,
    inventory_payload: bytes,
    timeout: float,
    show_progress: bool,
) -> tuple[Path, int]:
    cache_dir.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=f".{cache_dir.name}.staging-",
            dir=cache_dir.parent,
        )
    )
    progress = None
    try:
        if show_progress:
            from tqdm.auto import tqdm

            progress = tqdm(
                total=len(inventory.assets),
                desc="Cellucid web UI cache",
                unit="file",
            )
        downloaded_bytes = 0
        for asset in inventory.assets:
            response = _fetch_web_response(
                source_url=source_url,
                asset_path=asset.path,
                timeout=timeout,
            )
            _verify_fetched_asset(asset, response)
            target = stage / asset.path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(response.data)
            downloaded_bytes += len(response.data)
            if progress is not None:
                progress.update(1)
        (stage / WEB_ASSET_INVENTORY_FILENAME).write_bytes(inventory_payload)
        (stage / _WEB_CACHE_SOURCE_FILENAME).write_bytes(
            _source_control_bytes(source_url)
        )
        verify_web_ui_cache(stage, expected_source_url=source_url)
        return stage, downloaded_bytes
    except BaseException as error:
        try:
            _remove_web_cache_dir(stage)
        except BaseException as cleanup_error:
            raise cleanup_error from error
        raise
    finally:
        if progress is not None:
            progress.close()


def _activate_staged_generation(stage: Path, cache_dir: Path) -> None:
    backup = cache_dir.parent / f".{cache_dir.name}.previous-{uuid.uuid4().hex}"
    if not _require_web_cache_directory_or_missing(stage):
        raise FileNotFoundError(f"Staged web cache directory does not exist: {stage}")
    if _require_web_cache_directory_or_missing(backup):
        raise FileExistsError(f"Web cache backup path already exists: {backup}")
    had_active = _require_web_cache_directory_or_missing(cache_dir)
    if had_active:
        cache_dir.rename(backup)
    try:
        if _require_web_cache_directory_or_missing(cache_dir):
            raise FileExistsError(f"Web cache activation path already exists: {cache_dir}")
        if not _require_web_cache_directory_or_missing(stage):
            raise FileNotFoundError(f"Staged web cache directory does not exist: {stage}")
        stage.rename(cache_dir)
    except BaseException as activation_error:
        if had_active:
            try:
                if _require_web_cache_directory_or_missing(cache_dir):
                    raise FileExistsError(
                        f"Cannot restore the prior web cache over an existing path: {cache_dir}"
                    )
                if not _require_web_cache_directory_or_missing(backup):
                    raise FileNotFoundError(
                        f"Prior web cache directory is unavailable for restoration: {backup}"
                    )
                backup.rename(cache_dir)
            except BaseException as restoration_error:
                raise restoration_error from activation_error
        raise
    if had_active:
        _remove_web_cache_dir(backup)


def get_web_cache_dir() -> Path:
    """Return the exact current web cache directory."""
    return _web_cache_dir()


def clear_web_cache(*, cache_dir: str | Path | None = None) -> Path:
    """Delete the selected cache generation and propagate any cleanup failure."""
    resolved = _require_cache_dir(cache_dir)
    with _CACHE_LOCK:
        _require_web_cache_directory_or_missing(resolved)
        _remove_web_cache_dir(resolved)
    return resolved


def ensure_web_ui_cached(
    *,
    cache_dir: str | Path | None = None,
    source_url: str = CELLUCID_WEB_URL,
    force: bool = True,
    show_progress: bool = True,
    timeout: float = 15.0,
) -> WebCachePrefetchSummary:
    """Establish the source build, or verify only when ``force=False``."""
    if type(force) is not bool:
        raise TypeError("force must be a boolean")
    if type(show_progress) is not bool:
        raise TypeError("show_progress must be a boolean")
    source = _require_source_url(source_url)
    resolved = _require_cache_dir(cache_dir)
    if not isinstance(timeout, int | float) or isinstance(timeout, bool) or timeout <= 0:
        raise ValueError("timeout must be a positive number")

    if not force:
        with _CACHE_LOCK:
            cached_inventory = verify_web_ui_cache(
                resolved,
                expected_source_url=source,
            )
            return WebCachePrefetchSummary(
                cache_dir=resolved,
                source_url=source,
                build_id=cached_inventory.build_id,
                downloaded_files=0,
                downloaded_bytes=0,
                skipped_files=len(cached_inventory.assets),
            )

    inventory_response = _require_fetched_response(
        _fetch_web_response(
            source_url=source,
            asset_path=WEB_ASSET_INVENTORY_FILENAME,
            timeout=timeout,
        ),
        label="Web asset inventory",
    )
    if inventory_response.content_type != "application/json; charset=utf-8":
        raise ValueError(
            "Web asset inventory response Content-Type must be "
            "'application/json; charset=utf-8'"
        )
    inventory = parse_web_asset_inventory(inventory_response.data)

    with _CACHE_LOCK:
        _require_web_cache_directory_or_missing(resolved)
        stage, downloaded_bytes = _establish_staged_generation(
            cache_dir=resolved,
            source_url=source,
            inventory=inventory,
            inventory_payload=inventory_response.data,
            timeout=timeout,
            show_progress=show_progress,
        )
        try:
            _require_web_cache_directory_or_missing(resolved)
            _activate_staged_generation(stage, resolved)
        except BaseException:
            if stage.exists() or stage.is_symlink():
                _remove_web_cache_dir(stage)
            raise

    return WebCachePrefetchSummary(
        cache_dir=resolved,
        source_url=source,
        build_id=inventory.build_id,
        downloaded_files=len(inventory.assets),
        downloaded_bytes=downloaded_bytes,
        skipped_files=0,
    )
