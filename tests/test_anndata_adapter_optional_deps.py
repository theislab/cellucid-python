from __future__ import annotations

from pathlib import Path

import pytest


def test_anndata_adapter_zarr_missing_dependency_has_actionable_error(tmp_path: Path) -> None:
    """
    If `zarr` isn't installed, loading a `.zarr` store should fail with a clear,
    actionable error message (instead of a confusing stack trace).
    """
    import importlib.util

    if importlib.util.find_spec("zarr") is not None:
        pytest.skip("zarr is installed in this environment")

    # Minimal zarr "shape": presence of .zgroup is enough for adapter detection.
    zarr_dir = tmp_path / "demo.zarr"
    zarr_dir.mkdir()
    (zarr_dir / ".zgroup").write_text("{}", encoding="utf-8")

    from cellucid.anndata_adapter import AnnDataAdapter

    with pytest.raises(ImportError) as excinfo:
        AnnDataAdapter.from_file(zarr_dir)

    assert "pip install zarr" in str(excinfo.value)
