from __future__ import annotations

import hashlib
import json
import stat
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest


def _make_directory_symlink(link: Path, target: Path) -> None:
    try:
        link.symlink_to(target, target_is_directory=True)
    except (NotImplementedError, OSError) as error:
        pytest.skip(f"Directory symlinks are unavailable on this platform: {error}")


def _snapshot(directory: Path) -> dict[str, bytes]:
    return {
        path.relative_to(directory).as_posix(): path.read_bytes()
        for path in directory.rglob("*")
        if path.is_file()
    }


def _one_asset_responses() -> dict[str, tuple[bytes, str]]:
    index = (
        b"<!doctype html><html><head>"
        b'<meta name="cellucid-web-build-id" content="confinement-build">'
        b"</head><body></body></html>"
    )
    inventory = json.dumps(
        {
            "version": 1,
            "build_id": "confinement-build",
            "assets": [
                {
                    "path": "index.html",
                    "sha256": hashlib.sha256(index).hexdigest(),
                    "bytes": len(index),
                    "content_type": "text/html; charset=utf-8",
                }
            ],
        },
        separators=(",", ":"),
    ).encode()
    return {
        "cellucid-web-assets.json": (
            inventory,
            "application/json; charset=utf-8",
        ),
        "index.html": (index, "text/html; charset=utf-8"),
    }


@pytest.mark.parametrize("operation", ["clear", "force_refresh"])
@pytest.mark.parametrize("selection_mode", ["explicit", "default"])
@pytest.mark.parametrize(
    ("target_kind", "exception_type", "message"),
    [
        ("directory_symlink", ValueError, "symbolic link"),
        ("broken_symlink", ValueError, "symbolic link"),
        ("file", NotADirectoryError, "directory"),
    ],
)
def test_web_cache_mutations_reject_unsafe_selected_path_before_any_action(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    operation: str,
    selection_mode: str,
    target_kind: str,
    exception_type: type[Exception],
    message: str,
) -> None:
    from cellucid import web_cache

    selected = tmp_path / "selected-cache"
    target: Path | None = None
    expected_target_files: dict[str, bytes] | None = None
    expected_file_bytes: bytes | None = None

    if target_kind == "directory_symlink":
        target = tmp_path / "valuable-data"
        (target / "nested").mkdir(parents=True)
        (target / "keep.bin").write_bytes(b"\x00valuable\xff")
        (target / "nested" / "exact.txt").write_bytes(b"must survive exactly")
        expected_target_files = _snapshot(target)
        _make_directory_symlink(selected, target)
    elif target_kind == "broken_symlink":
        target = tmp_path / "missing-target"
        _make_directory_symlink(selected, target)
        assert not target.exists()
    elif target_kind == "file":
        expected_file_bytes = b"\x00this is not a cache directory\xff"
        selected.write_bytes(expected_file_bytes)
    else:
        raise AssertionError(f"Unhandled target kind: {target_kind}")

    if selection_mode == "default":
        monkeypatch.setattr(web_cache, "_web_cache_dir", lambda: selected)
        cache_arguments: dict[str, Path] = {}
    else:
        cache_arguments = {"cache_dir": selected}

    fetch = mock.Mock(side_effect=AssertionError("network must not be called"))
    monkeypatch.setattr(web_cache, "_fetch_web_response", fetch)

    with pytest.raises(exception_type, match=message) as error:
        if operation == "clear":
            web_cache.clear_web_cache(**cache_arguments)
        else:
            web_cache.ensure_web_ui_cached(
                **cache_arguments,
                source_url="https://viewer.example",
                force=True,
                show_progress=False,
            )

    assert str(selected) in str(error.value)
    fetch.assert_not_called()

    if target_kind == "directory_symlink":
        assert selected.is_symlink()
        assert target is not None
        assert expected_target_files is not None
        assert _snapshot(target) == expected_target_files
    elif target_kind == "broken_symlink":
        assert selected.is_symlink()
        assert target is not None
        assert not target.exists()
    else:
        assert selected.is_file()
        assert expected_file_bytes is not None
        assert selected.read_bytes() == expected_file_bytes


@pytest.mark.parametrize(
    ("mutation_fetch", "expected_calls"),
    [
        ("cellucid-web-assets.json", ["cellucid-web-assets.json"]),
        (
            "index.html",
            ["cellucid-web-assets.json", "index.html"],
        ),
    ],
)
def test_force_refresh_rechecks_selected_path_after_network_and_before_activation(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    mutation_fetch: str,
    expected_calls: list[str],
) -> None:
    from cellucid import web_cache

    selected = tmp_path / "selected-cache"
    target = tmp_path / "valuable-data"
    target.mkdir()
    (target / "keep.bin").write_bytes(b"\x00must survive refresh\xff")
    expected_target_files = _snapshot(target)
    responses = _one_asset_responses()
    calls: list[str] = []

    def fetch(*, source_url: str, asset_path: str, timeout: float):
        assert source_url == "https://viewer.example"
        assert timeout > 0
        calls.append(asset_path)
        if asset_path == mutation_fetch:
            _make_directory_symlink(selected, target)
        payload, content_type = responses[asset_path]
        return web_cache.FetchedWebResponse(
            data=payload,
            content_type=content_type,
            content_length=len(payload),
        )

    monkeypatch.setattr(web_cache, "_fetch_web_response", fetch)

    with pytest.raises(ValueError, match="symbolic link"):
        web_cache.ensure_web_ui_cached(
            cache_dir=selected,
            source_url="https://viewer.example",
            force=True,
            show_progress=False,
        )

    assert calls == expected_calls
    assert selected.is_symlink()
    assert _snapshot(target) == expected_target_files
    assert list(tmp_path.glob(".selected-cache.staging-*")) == []


def test_activation_rejects_selected_symlink_without_moving_either_directory(
    tmp_path: Path,
) -> None:
    from cellucid.web_cache import _activate_staged_generation

    target = tmp_path / "valuable-data"
    target.mkdir()
    (target / "keep.bin").write_bytes(b"\x00must survive activation\xff")
    expected_target_files = _snapshot(target)
    selected = tmp_path / "selected-cache"
    _make_directory_symlink(selected, target)
    stage = tmp_path / "complete-stage"
    stage.mkdir()
    (stage / "new.bin").write_bytes(b"new cache generation")
    expected_stage_files = _snapshot(stage)

    with pytest.raises(ValueError, match="symbolic link"):
        _activate_staged_generation(stage, selected)

    assert selected.is_symlink()
    assert _snapshot(target) == expected_target_files
    assert _snapshot(stage) == expected_stage_files


def test_windows_reparse_point_is_rejected_as_filesystem_indirection(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    from cellucid._server_base import _require_web_cache_directory_or_missing

    selected = tmp_path / "selected-cache"
    metadata = SimpleNamespace(
        st_mode=stat.S_IFDIR,
        st_file_attributes=stat.FILE_ATTRIBUTE_REPARSE_POINT,
    )
    monkeypatch.setattr(Path, "lstat", lambda _path: metadata)

    with pytest.raises(ValueError, match="reparse point"):
        _require_web_cache_directory_or_missing(selected)


def test_clear_web_cache_preserves_lexical_path_through_symlinked_parent(
    tmp_path: Path,
) -> None:
    from cellucid.web_cache import clear_web_cache

    real_parent = tmp_path / "real-parent"
    real_parent.mkdir()
    parent_alias = tmp_path / "parent-alias"
    _make_directory_symlink(parent_alias, real_parent)
    selected = parent_alias / "selected-cache"
    selected.mkdir()
    (selected / "cache.bin").write_bytes(b"selected cache bytes")

    cleared = clear_web_cache(cache_dir=selected)

    assert cleared == selected
    assert parent_alias.is_symlink()
    assert not (real_parent / "selected-cache").exists()
