"""Contract for the directories the web cache refuses to claim as ``cache_dir``.

``clear_web_cache`` removes the selected directory outright, and
``ensure_web_ui_cached(force=True)`` renames the whole existing directory aside
and then removes it, so both destroy every file the named directory holds. A
directory nobody set aside for the web build -- the filesystem root, the
current working directory, the home directory, or the directory holding every
home -- must therefore be refused before any path is inspected, fetched,
staged, renamed, or removed, whether it is named directly, through ``~``, or
through a symbolic link anywhere in the path.

The guard sits in ``_require_cache_dir``, the one function every entry point
that takes a ``cache_dir`` passes through, and it refuses exactly the set
``prepare()`` refuses as ``out_dir`` (see
``test_protected_output_directory_contract.py``) -- the two are one rule and
must not drift apart.
"""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path

import pytest

from cellucid import web_cache

USER_FILE_NAME = "USER_FILE.txt"
USER_FILE_TEXT = "work the user never set aside for a web cache\n"
SOURCE_URL = "https://viewer.example"
BUILD_ID = "protected-cache-build"
REFUSAL = "cache_dir must name a dedicated web cache directory"


@contextmanager
def _working_directory(directory: Path) -> Iterator[None]:
    previous = Path.cwd()
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(previous)


def _populate(directory: Path) -> None:
    """Fill a directory with the user work a cache replacement destroys."""
    directory.mkdir(parents=True, exist_ok=True)
    (directory / USER_FILE_NAME).write_text(USER_FILE_TEXT, encoding="utf-8")
    nested = directory / "my_notes"
    nested.mkdir(exist_ok=True)
    (nested / "note.md").write_text("notes\n", encoding="utf-8")


def _tree(directory: Path) -> set[str]:
    return {path.relative_to(directory).as_posix() for path in directory.rglob("*")}


def _fake_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point ``~`` at a scratch directory so no real home is ever named."""
    home = tmp_path / "homes" / "someone"
    home.mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    monkeypatch.delenv("HOMEDRIVE", raising=False)
    monkeypatch.delenv("HOMEPATH", raising=False)
    return home


def _make_directory_symlink(link: Path, target: Path) -> None:
    try:
        link.symlink_to(target, target_is_directory=True)
    except (NotImplementedError, OSError) as error:
        pytest.skip(f"Directory symlinks are unavailable on this platform: {error}")


def _published_build() -> dict[str, tuple[bytes, str]]:
    """One complete, verifiable web build served without touching the network."""
    index = (
        "<!doctype html><html><head>"
        f'<meta name="cellucid-web-build-id" content="{BUILD_ID}">'
        "</head><body></body></html>"
    ).encode()
    inventory = json.dumps(
        {
            "version": 1,
            "build_id": BUILD_ID,
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
        web_cache.WEB_ASSET_INVENTORY_FILENAME: (
            inventory,
            "application/json; charset=utf-8",
        ),
        "index.html": (index, "text/html; charset=utf-8"),
    }


@pytest.fixture
def source(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Serve a real build offline and record every asset the code asks for."""
    responses = _published_build()
    calls: list[str] = []

    def fetch(
        *,
        source_url: str,
        asset_path: str,
        timeout: float,
    ) -> web_cache.FetchedWebResponse:
        assert source_url == SOURCE_URL
        calls.append(asset_path)
        payload, content_type = responses[asset_path]
        return web_cache.FetchedWebResponse(
            data=payload,
            content_type=content_type,
            content_length=len(payload),
        )

    monkeypatch.setattr(web_cache, "_fetch_web_response", fetch)
    return calls


def _clear(cache_dir: object) -> object:
    return web_cache.clear_web_cache(cache_dir=cache_dir)  # type: ignore[arg-type]


def _force_refresh(cache_dir: object) -> object:
    return web_cache.ensure_web_ui_cached(
        cache_dir=cache_dir,  # type: ignore[arg-type]
        source_url=SOURCE_URL,
        force=True,
        show_progress=False,
    )


Sink = Callable[[object], object]

destructive_sink = pytest.mark.parametrize(
    "sink",
    [_clear, _force_refresh],
    ids=["clear_web_cache", "ensure_web_ui_cached(force=True)"],
)


# ---------------------------------------------------------------------------
# Refused targets
# ---------------------------------------------------------------------------


@destructive_sink
def test_the_working_directory_named_by_absolute_path_is_refused(
    tmp_path: Path,
    sink: Sink,
    source: list[str],
) -> None:
    """The reported data loss: ``cache_dir=os.getcwd()`` erased the user's files."""
    workdir = tmp_path / "workdir"
    _populate(workdir)
    before = _tree(workdir)

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        sink(os.getcwd())

    assert REFUSAL in str(error.value)
    assert _tree(workdir) == before
    assert (workdir / USER_FILE_NAME).read_text(encoding="utf-8") == USER_FILE_TEXT
    assert source == []


@destructive_sink
def test_the_home_directory_named_by_tilde_is_refused(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sink: Sink,
    source: list[str],
) -> None:
    """``~`` has a non-empty leaf, so the old root-only guard let it through."""
    home = _fake_home(tmp_path, monkeypatch)
    _populate(home)
    before = _tree(home)
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        sink("~")

    assert REFUSAL in str(error.value)
    assert _tree(home) == before
    assert (home / USER_FILE_NAME).read_text(encoding="utf-8") == USER_FILE_TEXT
    assert source == []


@destructive_sink
def test_the_home_directory_named_by_absolute_path_is_refused(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sink: Sink,
    source: list[str],
) -> None:
    home = _fake_home(tmp_path, monkeypatch)
    _populate(home)
    before = _tree(home)
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        sink(home)

    assert REFUSAL in str(error.value)
    assert _tree(home) == before
    assert source == []


@destructive_sink
def test_the_directory_that_holds_every_home_is_refused(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sink: Sink,
    source: list[str],
) -> None:
    """Replacing ``/Users`` or ``/home`` would remove every user's home."""
    home = _fake_home(tmp_path, monkeypatch)
    home_parent = home.parent
    _populate(home_parent / "another_user")
    before = _tree(home_parent)
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        sink(home_parent)

    assert REFUSAL in str(error.value)
    assert _tree(home_parent) == before
    assert source == []


@destructive_sink
def test_the_filesystem_root_is_refused(
    tmp_path: Path,
    sink: Sink,
    source: list[str],
) -> None:
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        sink(Path(Path.cwd().anchor))

    assert REFUSAL in str(error.value)
    assert source == []


@destructive_sink
@pytest.mark.parametrize("relative", [".", "./", "..", "../"])
def test_relative_names_for_a_protected_directory_are_refused(
    tmp_path: Path,
    relative: str,
    sink: Sink,
    source: list[str],
) -> None:
    """``.`` is the working directory and ``..`` is the directory holding it."""
    parent = tmp_path / "parent"
    workdir = parent / "workdir"
    _populate(workdir)
    _populate(parent)
    before_parent = _tree(parent)

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        sink(relative)

    assert REFUSAL in str(error.value)
    assert _tree(parent) == before_parent
    assert source == []


@destructive_sink
def test_a_symbolic_link_to_the_home_directory_is_refused_as_the_home_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sink: Sink,
    source: list[str],
) -> None:
    """A protected directory reached through an alias is the same directory.

    The leaf-symlink guard already refused this shape, but for the wrong
    reason. The refusal must name the protected directory, so a caller who
    dereferences the link themselves does not simply retry and lose the data.
    """
    home = _fake_home(tmp_path, monkeypatch)
    _populate(home)
    before = _tree(home)
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    link = workdir / "home_link"
    _make_directory_symlink(link, home)

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        sink(link)

    assert REFUSAL in str(error.value)
    assert str(home.resolve()) in str(error.value)
    assert link.is_symlink()
    assert _tree(home) == before
    assert (home / USER_FILE_NAME).read_text(encoding="utf-8") == USER_FILE_TEXT
    assert source == []


@destructive_sink
def test_a_symbolic_link_in_the_parent_chain_is_refused(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    sink: Sink,
    source: list[str],
) -> None:
    """The abspath-versus-resolve hole: only the leaf link used to be seen.

    ``<alias>/<user>`` has a real directory as its leaf, so the symbolic-link
    guard passes it, and a lexical absolute path never compares equal to the
    home directory it actually names.
    """
    home = _fake_home(tmp_path, monkeypatch)
    _populate(home)
    before = _tree(home)
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    alias = workdir / "homes_alias"
    _make_directory_symlink(alias, home.parent)
    through_alias = alias / home.name
    assert not through_alias.is_symlink()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        sink(through_alias)

    assert REFUSAL in str(error.value)
    assert _tree(home) == before
    assert (home / USER_FILE_NAME).read_text(encoding="utf-8") == USER_FILE_TEXT
    assert source == []


def test_every_cache_dir_entry_point_refuses_the_same_directory(
    tmp_path: Path,
    source: list[str],
) -> None:
    """One guard in ``_require_cache_dir`` covers every route to the sinks."""
    workdir = tmp_path / "workdir"
    _populate(workdir)
    before = _tree(workdir)

    with _working_directory(workdir):
        selected = os.getcwd()
        for call in (
            lambda: web_cache.clear_web_cache(cache_dir=selected),
            lambda: web_cache.ensure_web_ui_cached(
                cache_dir=selected,
                source_url=SOURCE_URL,
                force=True,
                show_progress=False,
            ),
            lambda: web_cache.ensure_web_ui_cached(
                cache_dir=selected,
                source_url=SOURCE_URL,
                force=False,
                show_progress=False,
            ),
            lambda: web_cache.verify_web_ui_cache(selected),
            lambda: web_cache.read_cached_web_asset(selected, "index.html"),
        ):
            with pytest.raises(ValueError) as error:
                call()
            assert REFUSAL in str(error.value)

    assert _tree(workdir) == before
    assert source == []


# ---------------------------------------------------------------------------
# The refusal itself
# ---------------------------------------------------------------------------


def test_refusal_names_the_argument_the_offending_path_and_the_fix(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        web_cache.clear_web_cache(cache_dir=os.getcwd())

    message = str(error.value)
    assert message.startswith("cache_dir must name a dedicated web cache directory, not ")
    assert str(workdir.resolve()) in message
    assert "removes the whole directory" in message
    assert "./cellucid-web-cache" in message


@destructive_sink
def test_refusal_creates_and_removes_nothing_beside_the_refused_directory(
    tmp_path: Path,
    sink: Sink,
    source: list[str],
) -> None:
    """The guard must precede the staging directory and the rename-aside backup."""
    parent = tmp_path / "parent"
    workdir = parent / "workdir"
    _populate(workdir)
    _populate(parent / "sibling")
    before_parent = _tree(parent)
    before_workdir = _tree(workdir)

    with _working_directory(workdir), pytest.raises(ValueError):
        sink(os.getcwd())

    assert _tree(parent) == before_parent
    assert _tree(workdir) == before_workdir
    assert list(parent.glob(".workdir.staging-*")) == []
    assert list(parent.glob(".workdir.previous-*")) == []
    assert source == []


def test_a_cache_dir_that_is_not_a_path_is_reported_as_such(tmp_path: Path) -> None:
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(TypeError) as error:
        web_cache.clear_web_cache(cache_dir=1)  # type: ignore[arg-type]

    assert str(error.value) == "cache_dir must be a native string or os.PathLike path."


# ---------------------------------------------------------------------------
# Accepted targets
# ---------------------------------------------------------------------------


def test_a_dedicated_cache_directory_still_publishes_verifies_reads_and_clears(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    source: list[str],
) -> None:
    """The guard must not cost the ordinary cache lifecycle."""
    home = _fake_home(tmp_path, monkeypatch)
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir):
        for selected in (
            workdir / "cellucid-web-cache",
            "./relative-cache",
            "~/tilde-cache",
        ):
            summary = web_cache.ensure_web_ui_cached(
                cache_dir=selected,
                source_url=SOURCE_URL,
                force=True,
                show_progress=False,
            )
            assert summary.build_id == BUILD_ID
            assert summary.downloaded_files == 1

            assert web_cache.verify_web_ui_cache(selected).build_id == BUILD_ID
            read = web_cache.read_cached_web_asset(selected, "index.html")
            assert read is not None
            assert BUILD_ID.encode() in read[0]

            cleared = web_cache.clear_web_cache(cache_dir=selected)
            assert not cleared.exists()

    assert (workdir / "relative-cache").exists() is False
    assert home.is_dir()


def test_cache_dir_none_still_resolves_to_the_default_location(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    source: list[str],
) -> None:
    """Omitting ``cache_dir`` must still select the process-independent default."""
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir):
        # The real default is only resolved, never removed, by this assertion.
        assert web_cache._require_cache_dir(None) == Path(
            os.path.abspath(web_cache._web_cache_dir())
        )

        default = tmp_path / "default-cache"
        monkeypatch.setattr(web_cache, "_web_cache_dir", lambda: default)
        assert web_cache.get_web_cache_dir() == default

        summary = web_cache.ensure_web_ui_cached(
            source_url=SOURCE_URL,
            force=True,
            show_progress=False,
        )
        assert summary.cache_dir == default
        assert (default / "index.html").is_file()

        assert web_cache.clear_web_cache() == default
        assert not default.exists()
