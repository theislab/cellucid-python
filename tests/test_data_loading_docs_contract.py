from __future__ import annotations

import re
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


# A path that exists only on the machine the documentation was written on. The
# placeholders the guides teach with — ``/Users/you/…``, ``/home/you``,
# ``/Users/<username>/`` — are deliberately excluded: they are instructions, not
# captured output. Everything else naming a real home or scratch directory is a
# leak, and one shipped in the same commit that added the screenshot checklist.
_CAPTURE_HOST_PATH = re.compile(
    r"/(?:Users|home)/(?!you\b|<username>|<user>)[A-Za-z0-9._-]+/"
    r"|/private/tmp/"
    r"|\.pyenv/versions/"
)


def test_no_documentation_page_names_a_path_from_the_machine_it_was_written_on() -> None:
    """Committed docs must not carry the author's home, scratch, or pyenv paths.

    Screenshots hide this class of leak in pixels, where no scan reaches, which
    is why terminal transcripts ship as ``text`` blocks. Prose has no such
    excuse, so the same rule is enforced here mechanically. ``docs/**`` is
    walked directly rather than through Sphinx: the screenshot checklist is in
    ``exclude_patterns`` and so is never built, but it is still committed to a
    public repository.
    """
    offenders: list[str] = []
    for page in sorted(DOCS_ROOT.rglob("*.md")):
        if "_build" in page.parts:
            continue
        for number, line in enumerate(page.read_text(encoding="utf-8").splitlines(), 1):
            if _CAPTURE_HOST_PATH.search(line):
                offenders.append(f"{page.relative_to(REPOSITORY_ROOT)}:{number}")
    assert offenders == [], offenders


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
    assert "`generate_datasets.py`" in guide
    assert "one root-level" in normalized
    assert "there are no tests, validation harness" in normalized
    assert "workflows, `.github/` directory, or site assets" in normalized
    # The root enumeration went stale silently once already: it kept claiming
    # "that is all" after LICENSE and the policy files were added. Pin every
    # name it now lists so the next addition fails here instead of shipping.
    for policy_file in (
        "LICENSE",
        "CONTRIBUTING.md",
        "SECURITY.md",
        "SUPPORT.md",
        "CODE_OF_CONDUCT.md",
    ):
        assert f"`{policy_file}`" in guide
    assert "cross-references this guide" in normalized

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
    server_guide = _read(SERVER_GUIDE)

    # ``cellucid serve`` on a prepared export prints the *absolute* directory it
    # resolved (``server/_server.py`` calls ``print_detail("Path", str(self.data_dir))``),
    # so a photograph of this command necessarily bakes the capture host's path
    # into the image pixels, where no scan can find it. The transcripts therefore
    # ship as verbatim ``text`` blocks, never as screenshots.
    assert "```{figure} ../../../_static/screenshots/server/pancreas-cli-serve.png" not in (
        server_guide
    )
    assert "```{figure} ../../../_static/screenshots/server/serve-anndata-cli.png" not in (
        server_guide
    )
    assert "$ cellucid serve ./pancreas_export --no-browser" in server_guide
    assert "      Path: /path/to/pancreas_export" in server_guide
    assert "[1/3] Validating dataset..." in server_guide
    assert "3,696 cells" in server_guide
    assert "3,753 genes" in server_guide
    assert "486-file web cache" in server_guide
    assert "--web-source-url http://127.0.0.1:4183" in server_guide
    assert "/?source=remote" in server_guide
    assert "/?anndata=true" in server_guide

    # The direct-AnnData startup prints a different URL from the prepared one, and
    # both are transcribed, so neither can quietly become the other.
    assert "$ cellucid serve lung_atlas.h5ad \\" in server_guide
    assert "[1/4] Detecting format..." in server_guide
    assert "Mode: read-only backed h5ad" in server_guide

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
