"""The shipped ``py.typed`` marker must buy a downstream consumer real types.

``cellucid`` ships ``src/cellucid/py.typed`` and declares ``Typing :: Typed``.
Both are promises to a type checker in *someone else's* project, so every
assertion here is made from outside the package: a scratch consumer module is
type-checked in its own directory, exactly as a dependent project would be.

The package resolves its exports through a module-level ``__getattr__`` so that
``import cellucid`` stays free of numpy, pandas, scipy, and anndata. A static
checker cannot see through that, so the exports are additionally re-imported
under ``if TYPE_CHECKING:``. That block must never execute at runtime, which is
what ``test_importing_cellucid_stays_lazy`` pins.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

import cellucid

CONSUMER_MODULE = "uses_cellucid.py"
_REVEALED_RE = re.compile(
    rf'^{re.escape(CONSUMER_MODULE)}:(?P<line>\d+): note: Revealed type is "(?P<type>.*)"$'
)

# ``__version__`` is a plain module-level assignment, not a lazy export.
LAZY_EXPORTS = tuple(name for name in cellucid.__all__ if name != "__version__")

# Imported eagerly by ``cellucid`` itself, these would defeat the fast-CLI-startup
# rationale recorded in ``src/cellucid/__init__.py``.
HEAVY_DEPENDENCIES = (
    "numpy",
    "pandas",
    "scipy",
    "anndata",
    "zarr",
    "numcodecs",
    "h5py",
    "tqdm",
    "IPython",
)


def _run_mypy_on_consumer(source: str, tmp_path: Path) -> str:
    """Type-check ``source`` as a standalone project that depends on cellucid."""
    pytest.importorskip("mypy")
    consumer_dir = tmp_path / "consumer"
    consumer_dir.mkdir()
    (consumer_dir / CONSUMER_MODULE).write_text(source, encoding="utf-8")
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "mypy",
            "--cache-dir",
            str(tmp_path / "mypy_cache"),
            CONSUMER_MODULE,
        ],
        capture_output=True,
        text=True,
        cwd=consumer_dir,
        check=False,
    )
    return completed.stdout + completed.stderr


def test_every_public_export_has_a_real_type_for_a_consumer(tmp_path: Path) -> None:
    """No export may resolve to ``Any``; the lazy ``__getattr__`` used to erase them all."""
    lines = ["import cellucid", ""]
    revealed_at: dict[int, str] = {}
    for name in LAZY_EXPORTS:
        lines.append(f"reveal_type(cellucid.{name})  # noqa: F821")
        revealed_at[len(lines)] = name
    output = _run_mypy_on_consumer("\n".join(lines) + "\n", tmp_path)

    revealed: dict[str, str] = {}
    for line in output.splitlines():
        match = _REVEALED_RE.match(line.strip())
        if match is None:
            continue
        name = revealed_at.get(int(match.group("line")))
        if name is not None:
            revealed[name] = match.group("type")

    assert set(revealed) == set(LAZY_EXPORTS), (
        f"mypy did not reveal a type for every export; output:\n{output}"
    )
    erased = sorted(name for name, revealed_type in revealed.items() if revealed_type == "Any")
    assert erased == [], f"these exports are still Any for a consumer: {erased}"


def test_a_wrong_keyword_argument_is_a_consumer_visible_type_error(tmp_path: Path) -> None:
    """``prepare(nonexistent_argument=123)`` used to type-check clean."""
    output = _run_mypy_on_consumer(
        "import cellucid\n\ncellucid.prepare(nonexistent_argument=123)\n",
        tmp_path,
    )
    assert 'Unexpected keyword argument "nonexistent_argument" for "prepare"' in output
    assert "[call-arg]" in output


def test_importing_cellucid_stays_lazy() -> None:
    """The ``TYPE_CHECKING`` re-exports must not run, or import cost regresses."""
    probe = (
        "import json, sys\n"
        "import cellucid\n"
        "print(json.dumps({\n"
        f"    'heavy': sorted(m for m in {HEAVY_DEPENDENCIES!r} if m in sys.modules),\n"
        "    'submodules': sorted(m for m in sys.modules if m.startswith('cellucid.')),\n"
        "}))\n"
    )
    completed = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=True,
    )
    loaded = json.loads(completed.stdout)
    assert loaded["heavy"] == [], (
        f"import cellucid eagerly loaded heavy dependencies: {loaded['heavy']}"
    )
    assert loaded["submodules"] == [], (
        f"import cellucid eagerly loaded submodules: {loaded['submodules']}"
    )


@pytest.mark.parametrize("name", LAZY_EXPORTS)
def test_each_export_still_resolves_at_runtime_through_getattr(name: str) -> None:
    assert getattr(cellucid, name) is not None


def test_an_unknown_attribute_still_raises_attribute_error() -> None:
    with pytest.raises(AttributeError, match="has no attribute 'nonexistent_export'"):
        cellucid.nonexistent_export  # noqa: B018


def test_the_typed_marker_that_these_promises_rest_on_is_shipped() -> None:
    package_dir = Path(cellucid.__file__).parent
    assert (package_dir / "py.typed").is_file()
