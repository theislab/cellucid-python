"""Contract for user-supplied paths and rejected input types.

Three behaviours are pinned here:

1. ``out_dir`` / ``exports_dir`` default to ``./exports`` and are resolved
   against the working directory that is current when the function is called,
   never against the working directory that happened to be active at import
   time.
2. Every user-supplied path entry point expands a leading ``~``.
3. ``obs`` and ``var`` must be ``pandas.DataFrame`` objects, and a wrong type is
   reported as a ``TypeError`` naming the received type instead of being fed
   into a length comparison that produces a false row count.
"""

from __future__ import annotations

import inspect
import os
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import cellucid.prepare_data as prepare_data
from cellucid.jupyter import CellucidViewer
from cellucid.prepare_data import generate_datasets_manifest, prepare
from cellucid.server import CellucidServer, _list_exported_datasets

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
        "latent_space": embedding.copy(),
        "obs": pd.DataFrame({"score": [0.25, 0.5, 0.75]}),
        "X_umap_2d": embedding,
        "dataset_name": "Path contract",
        "dataset_id": "path-contract",
        "obs_categorical_dtype": "uint16",
        "centroid_min_points": 1,
    }
    kwargs.update(overrides)
    return kwargs


# ---------------------------------------------------------------------------
# 1. Lazy, relative default export directory
# ---------------------------------------------------------------------------


def test_export_directory_defaults_are_relative_and_not_frozen_at_import() -> None:
    """A frozen ``Path.cwd()`` default leaks the build machine path into docs."""
    assert not hasattr(prepare_data, "DEFAULT_EXPORT_DIR")

    for function, parameter in (
        (prepare, "out_dir"),
        (generate_datasets_manifest, "exports_dir"),
    ):
        default = inspect.signature(function).parameters[parameter].default
        assert isinstance(default, str), (function, default)
        assert not Path(default).is_absolute(), (function, default)
        assert Path(default) == Path("exports"), (function, default)


def test_prepare_default_out_dir_follows_the_call_time_working_directory(
    tmp_path: Path,
) -> None:
    """A notebook that chdirs must not keep writing to the old directory."""
    import_time_cwd = tmp_path / "import_time"
    call_time_cwd = tmp_path / "call_time"
    import_time_cwd.mkdir()
    call_time_cwd.mkdir()

    with _working_directory(import_time_cwd):
        # Simulate the import-time observation that used to be captured.
        pass

    with _working_directory(call_time_cwd):
        prepare(**_prepare_kwargs())

    assert (call_time_cwd / "exports" / "dataset_identity.json").is_file()
    assert not (import_time_cwd / "exports").exists()


def test_generate_datasets_manifest_default_follows_the_call_time_directory(
    tmp_path: Path,
) -> None:
    call_time_cwd = tmp_path / "call_time"
    call_time_cwd.mkdir()

    with _working_directory(call_time_cwd):
        prepare(**_prepare_kwargs(out_dir=Path("exports") / "path-contract"))
        manifest_path = generate_datasets_manifest(default_dataset="path-contract")
        resolved_manifest = manifest_path.resolve()

    assert not manifest_path.is_absolute()
    assert resolved_manifest == (call_time_cwd / "exports" / "datasets.json").resolve()
    assert resolved_manifest.is_file()


# ---------------------------------------------------------------------------
# 2. ``~`` expansion on every user-supplied path
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
        prepare(**_prepare_kwargs(out_dir="~/exports/path-contract"))

    assert (home / "exports" / "path-contract" / "dataset_identity.json").is_file()
    assert not (workdir / "~").exists()


def test_generate_datasets_manifest_expands_the_home_directory(
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
        prepare(**_prepare_kwargs(out_dir=home / "exports" / "path-contract"))
        manifest_path = generate_datasets_manifest(
            "~/exports",
            default_dataset="path-contract",
        )

    assert manifest_path == home / "exports" / "datasets.json"
    assert not (workdir / "~").exists()


def test_server_entry_points_expand_the_home_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    home = tmp_path / "home"
    home.mkdir()
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))

    dataset_dir = home / "exports" / "path-contract"
    prepare(**_prepare_kwargs(out_dir=dataset_dir))

    with _working_directory(workdir):
        assert [entry["id"] for entry in _list_exported_datasets("~/exports")] == [
            "path-contract"
        ]

        server = CellucidServer(
            "~/exports/path-contract",
            port=0,
            open_browser=False,
            quiet=True,
            serve_web_ui=False,
        )
        assert server.data_dir == dataset_dir.resolve()

    assert not (workdir / "~").exists()


def test_jupyter_viewer_expands_the_home_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    home = tmp_path / "home"
    home.mkdir()
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))

    missing = "~/exports/never-prepared"
    with _working_directory(workdir), pytest.raises(FileNotFoundError) as error:
        CellucidViewer(missing, auto_open=False)

    message = str(error.value)
    assert "~" not in message
    assert str(home / "exports" / "never-prepared") in message
    assert not (workdir / "~").exists()


def test_anndata_server_expands_the_home_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pytest.importorskip("anndata")
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
    package_root = Path(prepare_data.__file__).resolve().parent
    expected = {
        "prepare_data.py": ["Path(out_dir).expanduser()", "Path(exports_dir).expanduser()"],
        "server.py": [
            "Path(data_dir).expanduser()",
            "Path(data_dir).expanduser().resolve()",
            "Path(web_cache_dir).expanduser().resolve()",
        ],
        "jupyter.py": [
            "Path(data_dir).expanduser().resolve()",
            "Path(web_cache_dir).expanduser().resolve()",
        ],
        "anndata_server.py": [
            "Path(data).expanduser()",
            "Path(web_cache_dir).expanduser().resolve()",
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
        prepare(**_prepare_kwargs(obs=obs, out_dir=tmp_path / "generation"))

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
            **_prepare_kwargs(
                var=var,
                gene_expression=np.ones((N_CELLS, 2), dtype=np.float32),
                out_dir=tmp_path / "generation",
            )
        )

    assert str(error.value) == f"var must be a pandas DataFrame, got {type_name}."


def test_missing_obs_and_var_still_raise_the_required_value_errors(
    tmp_path: Path,
) -> None:
    kwargs = _prepare_kwargs(out_dir=tmp_path / "no-obs")
    kwargs["obs"] = None
    with pytest.raises(ValueError, match="obs DataFrame is required."):
        prepare(**kwargs)

    with pytest.raises(ValueError, match="var DataFrame must be provided"):
        prepare(
            **_prepare_kwargs(
                gene_expression=np.ones((N_CELLS, 2), dtype=np.float32),
                out_dir=tmp_path / "no-var",
            )
        )


def test_row_count_mismatch_still_reports_the_real_row_counts(tmp_path: Path) -> None:
    with pytest.raises(ValueError) as error:
        prepare(
            **_prepare_kwargs(
                obs=pd.DataFrame({"score": [0.1, 0.2]}),
                out_dir=tmp_path / "generation",
            )
        )

    assert str(error.value) == "obs has 2 rows, but embeddings have 3 cells."
