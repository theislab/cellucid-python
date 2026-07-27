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
        notebook = _load_notebook(path)
        text = json.dumps(notebook)
        assert "endocrinogenesis_day15.5_raw.h5ad" in text
        assert "/lustre/" not in text
        assert "he_developmental_complete" not in text
        assert "data/experiments/" not in text
        assert "bundle.cleanup" not in text


def test_jupyter_pancreas_executes_direct_anndata_diagnostics_before_cleanup() -> None:
    notebook = _load_notebook(TUTORIALS / "jupyter_pancreas.ipynb")
    source = "\n".join(
        "".join(cell["source"])
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )

    assert "show_anndata(" in source
    assert "assert adata.n_obs == 2531" in source
    assert 'viewer.set_color_by("clusters")' in source
    assert "viewer.debug_connection(timeout=15)" in source
    assert 'identity["stats"]["n_cells"]' in source
    assert 'identity["n_cells"]' not in source
    assert source.index("viewer.debug_connection") < source.index("viewer.stop()")
    assert "web_source_url" not in source
    assert "show(" not in source


def test_tutorial_notebook_inventory_contains_only_reproducible_pancreas_examples() -> None:
    notebook_names = sorted(path.name for path in TUTORIALS.glob("*.ipynb"))
    assert notebook_names == ["jupyter_pancreas.ipynb", "prepare_pancreas.ipynb"]

    recipe_gallery = (TUTORIALS / "04_real_world_dataset_recipes_gallery.md").read_text(
        encoding="utf-8"
    )
    jupyter_gallery = (
        TUTORIALS / "05_jupyter_embedding_hooks_sessions_gallery.md"
    ).read_text(encoding="utf-8")
    assert "prepare_pancreas" in recipe_gallery
    assert "jupyter_pancreas" in jupyter_gallery
    for stale_name in (
        "prepare_hdca",
        "prepare_he",
        "jupyter_hooks_sessions_he",
        "jupyter_session_bridge_he",
    ):
        assert stale_name not in recipe_gallery
        assert stale_name not in jupyter_gallery

    start_here = (TUTORIALS / "00_start_here.md").read_text(encoding="utf-8")
    success = start_here.split("## What success looks like", 1)[1].split(
        "## Troubleshooting", 1
    )[0]
    assert "{doc}`prepare_pancreas`" in success
    assert "{doc}`jupyter_pancreas`" in success
    assert "successful frontend round trip" in success


def test_vector_tutorial_uses_current_particle_flow_ui_language() -> None:
    tutorial = (
        TUTORIALS / "33_vector_fields_and_velocity_overlay_end_to_end.md"
    ).read_text(encoding="utf-8")
    lowered = tutorial.lower()
    assert "animated particle flow" in lowered
    assert "**visualization → vector field overlay**" in lowered
    assert "**show overlay**" in lowered
    assert "**particle density:**" in lowered
    assert "**flow speed:**" in lowered
    assert "**trail length:**" in lowered
    assert "overlay scale" not in lowered
    assert "arrow" not in lowered
