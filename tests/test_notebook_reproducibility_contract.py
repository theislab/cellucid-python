from __future__ import annotations

import json
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TUTORIALS = (
    REPOSITORY_ROOT
    / "docs"
    / "user_guide"
    / "python_package"
    / "f_notebooks_tutorials"
)
EXAMPLES = REPOSITORY_ROOT / "examples"
PANCREAS_H5AD = TUTORIALS / "datasets" / "endocrinogenesis_day15.5_raw.h5ad"
READER_NOTEBOOKS = (
    TUTORIALS / "jupyter_pancreas.ipynb",
    TUTORIALS / "prepare_pancreas.ipynb",
    EXAMPLES / "jupyter_cellucid_smoketest.ipynb",
)


def _load_notebook(path: Path) -> dict[str, object]:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(notebook, dict), path
    assert notebook.get("nbformat") == 4, path
    assert isinstance(notebook.get("metadata"), dict), path
    cells = notebook.get("cells")
    assert isinstance(cells, list) and cells, path

    cell_ids: list[str] = []
    for cell in cells:
        assert isinstance(cell, dict), path
        assert cell.get("cell_type") in {"code", "markdown", "raw"}, path
        assert isinstance(cell.get("metadata"), dict), path
        assert isinstance(cell.get("source"), list), path
        cell_id = cell.get("id")
        assert isinstance(cell_id, str) and cell_id, path
        cell_ids.append(cell_id)
        if cell.get("cell_type") == "code":
            assert cell.get("execution_count") is None, path
            assert cell.get("outputs") == [], path

    assert len(cell_ids) == len(set(cell_ids)), path
    return notebook


def test_reader_notebooks_are_valid_clean_json_with_public_checked_in_input() -> None:
    assert PANCREAS_H5AD.is_file()
    assert PANCREAS_H5AD.stat().st_size > 1_000_000

    for path in READER_NOTEBOOKS:
        _load_notebook(path)
