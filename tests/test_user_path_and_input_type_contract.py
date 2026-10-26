"""Contract for user-supplied paths and rejected input types.

Three behaviours are pinned here:

1. A relative ``out_path`` is resolved against the working directory that is
   current when :func:`cellucid.prepare` is called.
2. Every user-supplied path entry point expands a leading ``~``.
3. ``obs`` and ``var`` must be ``pandas.DataFrame`` objects, and a wrong type is
   reported as a ``TypeError`` naming the received type instead of being fed
   into a length comparison that produces a false row count.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cellucid
from cellucid import prepare

N_CELLS = 3


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
        "obs": pd.DataFrame({"score": [0.25, 0.5, 0.75]}),
        "embeddings": {"X_umap": embedding},
        "dataset_name": "Path contract",
        "dataset_id": "path-contract",
    }
    kwargs.update(overrides)
    return kwargs


# ---------------------------------------------------------------------------
# 1. A relative output path follows the call-time working directory
# ---------------------------------------------------------------------------


def test_prepare_resolves_a_relative_out_path_at_call_time(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    first.mkdir()
    second.mkdir()
    with _working_directory(second):
        written = prepare(
            "dataset.h5ad",
            format="h5ad",
            **_prepare_kwargs(
                gene_expression=np.ones((N_CELLS, 1), dtype=np.float32),
                var=pd.DataFrame(index=["gene-1"]),
            ),
        )

    assert written == (second / "dataset.h5ad").resolve()
    assert written.is_file()
    assert not (first / "dataset.h5ad").exists()


# ---------------------------------------------------------------------------
# 2. Home-directory expansion
# ---------------------------------------------------------------------------


def test_prepare_expands_the_home_directory_instead_of_creating_a_literal_tilde(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    home = tmp_path / "home"
    home.mkdir()
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))

    with _working_directory(workdir):
        written = prepare("~/prepared/path-contract.zarr", format="zarr", **_prepare_kwargs())

    assert written == (home / "prepared" / "path-contract.zarr").resolve()
    assert (written / ".zgroup").is_file()
    assert not (workdir / "~").exists()


def test_anndata_server_expands_the_home_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from cellucid.anndata_server import AnnDataServer

    home = tmp_path / "home"
    home.mkdir()
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))

    with _working_directory(workdir), pytest.raises(FileNotFoundError) as error:
        AnnDataServer(
            "~/never-prepared.h5ad",
            quiet=True,
            dataset_name="Path contract",
            dataset_id="path-contract",
        )

    message = str(error.value)
    assert "~" not in message
    assert str(home / "never-prepared.h5ad") in message


def test_every_user_supplied_path_entry_point_expands_the_home_directory() -> None:
    """Guard the whole class: the sites fixed here must keep calling expanduser."""
    package_root = Path(cellucid.__file__).resolve().parent
    expected = {
        "prepare_data/_prepare.py": [
            "Path(out_path).expanduser()",
        ],
        "jupyter/_viewer.py": [
            "Path(web_cache_dir).expanduser().resolve()",
        ],
        "anndata_server.py": [
            "Path(web_cache_dir).expanduser().resolve()",
        ],
        "_anndata_files.py": [
            "Path(path).expanduser().resolve()",
        ],
        "session_bundle.py": ["Path(path).expanduser().resolve()", "Path(dest).expanduser().resolve()"],
        "web_cache.py": ["cache_dir.expanduser()"],
    }
    for filename, fragments in expected.items():
        source = (package_root / filename).read_text(encoding="utf-8")
        for fragment in fragments:
            assert fragment in source, f"{filename} no longer expands: {fragment}"


# ---------------------------------------------------------------------------
# 3. obs / var type validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("obs", "type_name"),
    [
        ({"cell_type": ["a", "b", "c"]}, "dict"),
        (np.zeros((N_CELLS, 2), dtype=np.float32), "ndarray"),
        ([{"cell_type": "a"}], "list"),
    ],
)
def test_prepare_rejects_non_dataframe_obs_with_the_received_type(
    tmp_path: Path,
    obs: object,
    type_name: str,
) -> None:
    with pytest.raises(TypeError) as error:
        prepare(tmp_path / "generation.h5ad", format="h5ad", **_prepare_kwargs(obs=obs))

    assert str(error.value) == f"obs must be a pandas DataFrame, got {type_name}."


@pytest.mark.parametrize(
    ("var", "type_name"),
    [
        ({"gene": ["g1", "g2"]}, "dict"),
        (np.zeros((2, 1), dtype=np.float32), "ndarray"),
        (["g1", "g2"], "list"),
    ],
)
def test_prepare_rejects_non_dataframe_var_with_the_received_type(
    tmp_path: Path,
    var: object,
    type_name: str,
) -> None:
    with pytest.raises(TypeError) as error:
        prepare(
            tmp_path / "generation.h5ad",
            format="h5ad",
            **_prepare_kwargs(
                var=var,
                gene_expression=np.ones((N_CELLS, 2), dtype=np.float32),
            ),
        )

    assert str(error.value) == f"var must be a pandas DataFrame, got {type_name}."


def test_missing_gene_names_are_a_value_error_naming_the_remedy(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="The genes have no names"):
        prepare(
            tmp_path / "no-var.h5ad",
            format="h5ad",
            **_prepare_kwargs(gene_expression=np.ones((N_CELLS, 2), dtype=np.float32)),
        )
    assert list(tmp_path.iterdir()) == []


def test_row_count_mismatch_still_reports_the_real_row_counts(tmp_path: Path) -> None:
    with pytest.raises(ValueError) as error:
        prepare(
            tmp_path / "generation.h5ad",
            format="h5ad",
            **_prepare_kwargs(obs=pd.DataFrame({"score": [0.1, 0.2]})),
        )

    assert str(error.value) == "obs has 2 rows, but the embeddings have 3."
