"""Contract for the directories ``prepare()`` refuses to claim as ``out_dir``.

Publishing a generation renames the entire existing target aside and then
removes it, so ``force=True`` destroys every file the named directory holds.
A directory nobody set aside for one dataset -- the filesystem root, the
current working directory, the home directory, or the directory holding every
home -- must therefore be refused before any path is created, locked, written,
or removed, whether it is named directly, through ``~``, or through a symbolic
link.

The peer R package refuses the same set in ``.validate_output_path()``, so the
two writers agree on which directories are never an export target.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cellucid.prepare_data import prepare

N_CELLS = 3
USER_FILE_NAME = "USER_FILE.txt"
USER_FILE_TEXT = "work the user never set aside for an export\n"


@contextmanager
def _working_directory(directory: Path) -> Iterator[None]:
    previous = Path.cwd()
    os.chdir(directory)
    try:
        yield
    finally:
        os.chdir(previous)


def _prepare_kwargs(**overrides: object) -> dict[str, object]:
    embedding = np.array(
        [
            [-3.0, 1.0],
            [0.5, 5.0],
            [8.0, -2.0],
        ],
        dtype=np.float32,
    )
    kwargs: dict[str, object] = {
        "latent_space": embedding.copy(),
        "obs": pd.DataFrame({"score": [0.25, 0.5, 0.75]}),
        "X_umap_2d": embedding,
        "dataset_name": "Protected directory contract",
        "dataset_id": "protected-directory-contract",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
    }
    kwargs.update(overrides)
    return kwargs


def _populate(directory: Path) -> None:
    """Fill a directory with the user work a replacement generation destroys."""
    (directory / USER_FILE_NAME).write_text(USER_FILE_TEXT, encoding="utf-8")
    nested = directory / "my_notes"
    nested.mkdir()
    (nested / "note.md").write_text("notes\n", encoding="utf-8")


def _tree(directory: Path) -> set[str]:
    return {path.relative_to(directory).as_posix() for path in directory.rglob("*")}


def _fake_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point ``~`` at a scratch directory so no real home is ever named."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    monkeypatch.delenv("HOMEDRIVE", raising=False)
    monkeypatch.delenv("HOMEPATH", raising=False)
    return home


# ---------------------------------------------------------------------------
# Refused targets
# ---------------------------------------------------------------------------


def test_prepare_refuses_the_working_directory_named_by_absolute_path(
    tmp_path: Path,
) -> None:
    """The reported data loss: ``out_dir=os.getcwd()`` erased the user's files."""
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    _populate(workdir)
    before = _tree(workdir)

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        prepare(**_prepare_kwargs(out_dir=str(Path.cwd()), force=True))

    assert "out_dir must name a dedicated dataset output directory" in str(error.value)
    assert _tree(workdir) == before
    assert (workdir / USER_FILE_NAME).read_text(encoding="utf-8") == USER_FILE_TEXT


def test_prepare_refuses_the_home_directory_named_by_tilde(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """``~`` has a non-empty leaf, so a name-only guard lets it through."""
    home = _fake_home(tmp_path, monkeypatch)
    _populate(home)
    before = _tree(home)
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        prepare(**_prepare_kwargs(out_dir="~", force=True))

    assert "out_dir must name a dedicated dataset output directory" in str(error.value)
    assert _tree(home) == before
    assert (home / USER_FILE_NAME).read_text(encoding="utf-8") == USER_FILE_TEXT


def test_prepare_refuses_the_home_directory_named_by_absolute_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    home = _fake_home(tmp_path, monkeypatch)
    _populate(home)
    before = _tree(home)
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError):
        prepare(**_prepare_kwargs(out_dir=home, force=True))

    assert _tree(home) == before


def test_prepare_refuses_a_symbolic_link_to_the_home_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A protected directory reached through an alias is the same directory."""
    home = _fake_home(tmp_path, monkeypatch)
    _populate(home)
    before = _tree(home)
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    link = workdir / "home_link"
    link.symlink_to(home, target_is_directory=True)
    before_workdir = _tree(workdir)

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        prepare(**_prepare_kwargs(out_dir=link, force=True))

    # The refusal must come from the path guard, not from the symbolic-link
    # check inside the export lock, which has already created a lock file.
    assert "out_dir must name a dedicated dataset output directory" in str(error.value)
    assert _tree(workdir) == before_workdir
    assert link.is_symlink()
    assert _tree(home) == before
    assert (home / USER_FILE_NAME).read_text(encoding="utf-8") == USER_FILE_TEXT


def test_prepare_refuses_the_directory_that_holds_every_home(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Replacing ``/Users`` or ``/home`` would remove every user's home."""
    home = _fake_home(tmp_path, monkeypatch)
    home_parent = home.parent
    (home_parent / "another_user").mkdir()
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    before = _tree(home_parent)

    with _working_directory(workdir), pytest.raises(ValueError):
        prepare(**_prepare_kwargs(out_dir=home_parent, force=True))

    assert _tree(home_parent) == before


def test_prepare_refuses_the_filesystem_root(tmp_path: Path) -> None:
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        prepare(**_prepare_kwargs(out_dir=Path(Path.cwd().anchor), force=True))

    assert "out_dir must name a dedicated dataset output directory" in str(error.value)


@pytest.mark.parametrize("relative", [".", "..", "./", "../"])
def test_prepare_refuses_relative_names_for_a_protected_directory(
    tmp_path: Path,
    relative: str,
) -> None:
    """``.`` is the working directory and ``..`` is the directory holding it."""
    parent = tmp_path / "parent"
    workdir = parent / "workdir"
    workdir.mkdir(parents=True)
    _populate(workdir)
    _populate(parent)
    before_parent = _tree(parent)

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        prepare(**_prepare_kwargs(out_dir=relative, force=True))

    assert "out_dir must name a dedicated dataset output directory" in str(error.value)
    assert _tree(parent) == before_parent


def test_refusal_names_the_argument_the_offending_path_and_the_fix(
    tmp_path: Path,
) -> None:
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(ValueError) as error:
        prepare(**_prepare_kwargs(out_dir=str(Path.cwd())))

    message = str(error.value)
    assert message.startswith("out_dir must name a dedicated dataset output directory, not ")
    assert str(workdir.resolve()) in message
    assert "./exports/my_dataset" in message


def test_refusal_creates_nothing_beside_the_refused_directory(
    tmp_path: Path,
) -> None:
    """The guard must precede the parent ``mkdir`` and the sibling lock file."""
    parent = tmp_path / "parent"
    workdir = parent / "workdir"
    workdir.mkdir(parents=True)
    before_parent = _tree(parent)
    before_workdir = _tree(workdir)

    with _working_directory(workdir), pytest.raises(ValueError):
        prepare(**_prepare_kwargs(out_dir=str(Path.cwd()), force=True))

    assert _tree(parent) == before_parent
    assert _tree(workdir) == before_workdir


def test_prepare_rejects_an_out_dir_that_is_not_a_path(tmp_path: Path) -> None:
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    with _working_directory(workdir), pytest.raises(TypeError) as error:
        prepare(**_prepare_kwargs(out_dir=1))

    assert str(error.value) == "out_dir must be a native string or os.PathLike path."


# ---------------------------------------------------------------------------
# Accepted targets
# ---------------------------------------------------------------------------


def test_a_dedicated_child_directory_still_publishes_with_and_without_force(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The guard must not cost the ordinary export or the replacement path."""
    home = _fake_home(tmp_path, monkeypatch)
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    dataset_dir = workdir / "exports" / "protected-directory-contract"

    with _working_directory(workdir):
        prepare(**_prepare_kwargs(out_dir=dataset_dir))
        assert (dataset_dir / "dataset_identity.json").is_file()

        prepare(**_prepare_kwargs(out_dir=dataset_dir, force=True))
        assert (dataset_dir / "dataset_identity.json").is_file()

        prepare(**_prepare_kwargs(out_dir="~/exports/tilde-child"))
        assert (home / "exports" / "tilde-child" / "dataset_identity.json").is_file()

        prepare(**_prepare_kwargs(out_dir="./exports/relative-child"))
        assert (
            workdir / "exports" / "relative-child" / "dataset_identity.json"
        ).is_file()
