from __future__ import annotations

import struct
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPOSITORY_ROOT / "docs"
DATA_LOADING = DOCS_ROOT / "user_guide" / "web_app" / "b_data_loading"
CUSTOM_REPOSITORY = DATA_LOADING / "11_custom_dataset_repository.md"
SERVER_GUIDE = DATA_LOADING / "04_server_tutorial.md"
JUPYTER_GUIDE = DATA_LOADING / "05_jupyter_tutorial.md"
REFERENCE_VALUE = "theislab/cellucid-demo-custom-datasets/exports"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _png_dimensions(path: Path) -> tuple[int, int]:
    payload = path.read_bytes()
    assert payload[:8] == b"\x89PNG\r\n\x1a\n", path
    assert payload[12:16] == b"IHDR", path
    return struct.unpack(">II", payload[16:24])


def test_central_navigation_starts_with_only_top_level_sections_open() -> None:
    source = _read(DOCS_ROOT / "conf.py")
    assert source.count('"show_nav_level": 1,') == 1
    assert '"show_nav_level": 2,' not in source
    assert 'html_theme = "pydata_sphinx_theme"' in source


def test_custom_repository_guide_publishes_the_exact_current_catalog_contract() -> None:
    guide = _read(CUSTOM_REPOSITORY)
    normalized = " ".join(guide.split())

    assert f"```text\n   {REFERENCE_VALUE}\n   ```" in guide
    assert "resolves deterministically to the `main` branch" in guide
    assert "does not probe `master`, `gh-pages`, or any other branch" in guide
    assert "owner/repo@release-2026/path/to/exports" in guide
    assert "owner/repo/release-2026/path/to/exports" in guide

    for dataset_id, n_cells, n_genes in (
        ("synthetic-cell-types-2d", 72, 8),
        ("synthetic-development-3d", 96, 10),
        ("synthetic-trajectory-1d", 48, 6),
    ):
        assert f"`{dataset_id}`" in guide
        assert f"{n_cells} cells, {n_genes} genes" in guide
        assert f"?github={REFERENCE_VALUE}&dataset={dataset_id}" in guide

    for exact_catalog_rule in (
        "`version`, `default`, and `datasets`",
        "`version` is exactly `1`",
        "unique `id`",
        "safe relative",
        "ends in `/`",
        "`name`, `description`, `n_cells`, and `n_genes`",
        "exactly matches the canonical value",
        "validates `datasets.json` and every",
        "listed `dataset_identity.json`",
    ):
        assert exact_catalog_rule in normalized

    assert 'generate_datasets_manifest(\n    "./exports",\n    default_dataset="study-a",' in guide
    assert "cellucid serve ./exports --no-browser" in guide
    assert "/?source=remote" in guide
    assert (
        "https://raw.githubusercontent.com/theislab/"
        "cellucid-demo-custom-datasets/main/exports/datasets.json"
    ) in guide


def test_loading_guides_do_not_reintroduce_stale_branch_or_server_routes() -> None:
    combined = "\n".join(_read(path) for path in sorted(DATA_LOADING.glob("*.md")))
    for stale_text in (
        "try common defaults",
        "owner/repo/gh-pages/exports",
        "http://127.0.0.1:8765/health",
        'generate_datasets_manifest("./exports")',
        "Dataset Connections",
        "Browse local data",
    ):
        assert stale_text not in combined

    assert "http://127.0.0.1:8765/_cellucid/health" in combined
    assert "http://127.0.0.1:8765/_cellucid/datasets" in combined
    assert "**Session**" in combined
    assert "**Local data:**" in combined
    assert "**Prepared**" in combined
    assert "**H5AD**" in combined
    assert "**Zarr ZIP**" in combined


def test_real_pancreas_server_and_jupyter_evidence_is_wired_into_the_guides() -> None:
    screenshot = DOCS_ROOT / "_static" / "screenshots" / "server" / "pancreas-cli-serve.png"
    assert _png_dimensions(screenshot) == (1440, 1000)

    server_guide = _read(SERVER_GUIDE)
    assert "```{figure} ../../../_static/screenshots/server/pancreas-cli-serve.png" in server_guide
    assert ":width: 1440px" in server_guide
    assert "3,696 cells" in server_guide
    assert "3,753 genes" in server_guide
    assert "429-file web cache" in server_guide
    assert "--web-source-url http://127.0.0.1:4173" in server_guide
    assert "/?source=remote" in server_guide
    assert "/?anndata=true" in server_guide

    jupyter_guide = _read(JUPYTER_GUIDE)
    assert "Pancreas Jupyter walkthrough" in jupyter_guide
    assert (
        "../../python_package/f_notebooks_tutorials/05_jupyter_embedding_hooks_sessions_gallery"
    ) in jupyter_guide

    jupyter_gallery = _read(
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "f_notebooks_tutorials"
        / "05_jupyter_embedding_hooks_sessions_gallery.md"
    )
    assert "jupyter_pancreas" in jupyter_gallery
    assert "pancreas-notebook-embed.png" in jupyter_gallery
    assert "pancreas-debug-connection.png" in jupyter_gallery
    assert "show_anndata(DATA_PATH" in jupyter_gallery
    assert "normal released-build resolution path" in jupyter_gallery
    assert "show_anndata(adata)" not in jupyter_gallery
    assert "no-argument" not in jupyter_gallery

    overview = _read(DATA_LOADING / "01_loading_options_overview.md")
    local_demo = _read(DATA_LOADING / "02_local_demo_tutorial.md")
    assert "Browser **Prepared** picker" in overview
    assert "**Prepared** file picker" in local_demo
    assert "Browser **Folder** picker" not in overview
    assert "**Folder** file picker" not in local_demo
    assert "I clicked Folder" not in _read(
        DATA_LOADING / "03_browser_file_picker_tutorial.md"
    )
    assert "deliberately exports only the checked-in source file's real 2-D UMAP" in local_demo
    assert "shows how the real Pancreas dataset creates the 1-D" not in local_demo


def test_python_and_r_sharing_navigation_points_to_the_reference_repository() -> None:
    paths = (
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "d_viewing_apis"
        / "03_choose_your_workflow_decision_tree.md",
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "d_viewing_apis"
        / "07_exported_directory_mode_show_and_serve.md",
        DOCS_ROOT
        / "user_guide"
        / "r_package"
        / "d_viewing_loading"
        / "02_host_exports_for_sharing.md",
    )
    for path in paths:
        assert "11_custom_dataset_repository" in _read(path), path

    r_guide = _read(paths[-1])
    assert REFERENCE_VALUE in r_guide
    assert "`%||%`" not in r_guide
    assert "if (!file.exists(ident_path)) next" not in r_guide
