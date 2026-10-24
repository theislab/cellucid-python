from __future__ import annotations

import re
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS_ROOT = REPOSITORY_ROOT / "docs"
DATA_LOADING = DOCS_ROOT / "user_guide" / "web_app" / "b_data_loading"
CUSTOM_REPOSITORY = DATA_LOADING / "11_custom_dataset_repository.md"
SERVER_GUIDE = DATA_LOADING / "04_server_tutorial.md"
JUPYTER_GUIDE = DATA_LOADING / "05_jupyter_tutorial.md"
REFERENCE_VALUE = "theislab/cellucid-demo-custom-datasets/catalog"


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

    # Shorthand reads the default branch as GitHub resolves HEAD: no probing.
    assert "reads the repository's default branch (`main`, `master` or any other)" in normalized
    assert "Cellucid never probes branch names" in normalized
    assert "{doc}`15_host_on_github`" in guide
    assert "owner/repo@release-2026/path/to/catalog" in guide
    assert "owner/repo/release-2026/path/to/catalog" in guide
    assert "`generate_datasets.py`" in guide

    for dataset_id, n_cells, n_genes in (
        ("synthetic-cell-types-2d", 72, 8),
        ("synthetic-development-3d", 96, 10),
        ("synthetic-trajectory-1d", 48, 6),
    ):
        assert f"`{dataset_id}`" in guide
        assert f"{n_cells} cells, {n_genes} genes" in guide

    # The catalog lists AnnData files (schema 2); names, counts and citations
    # come from each file, never from the catalog.
    assert '"schema": 2' in guide
    for exact_catalog_rule in (
        "`schema`, `default`, and `datasets`",
        "`schema` is exactly `2`",
        "Every `id` is unique and follows the portable id rule",
        "names a Zarr store folder (ending in `/`), an `.h5ad`, or a `.zarr.zip`",
        "with its exact `bytes` and lowercase `sha256`",
        "Names, descriptions, citations and counts are not in the catalog",
    ):
        assert exact_catalog_rule in normalized
    for removed in (
        "`version`",
        "dataset_identity.json",
        "generate_datasets_manifest",
        "/?source=remote",
        "**Prepared**",
    ):
        assert removed not in guide


def test_loading_guides_do_not_reintroduce_stale_branch_or_server_routes() -> None:
    combined = "\n".join(_read(path) for path in sorted(DATA_LOADING.glob("*.md")))
    for stale_text in (
        "try common defaults",
        "owner/repo/gh-pages/exports",
        "http://127.0.0.1:8765/health",
        "generate_datasets_manifest",
        "Dataset Connections",
        "Browse local data",
        "**Prepared**",
        "exportsBaseUrl",
        "state-snapshots.json",
        "/?source=remote",
    ):
        assert stale_text not in combined

    assert "/_cellucid/health" in combined
    assert "/_cellucid/datasets" in combined
    assert "**Session**" in combined
    assert "**Open…**" in combined
    assert "**File…**" in combined
    assert "**Zarr folder…**" in combined
    assert "**Link or server**" in combined
    # The Session panel's earlier controls are gone: the Files buttons, the
    # Server and Open from link fields with Connect, and the Transfers list.
    for removed in (
        "**Files**",
        "**H5AD**",
        "**Zarr ZIP**",
        "**Zarr folder**",
        "Open from link**",
        "**Server** field",
        "Reconnect",
        "**Transfers**",
    ):
        assert removed not in combined, removed
    assert "?sampleCatalog=" in combined


def test_real_pancreas_server_and_jupyter_evidence_is_wired_into_the_guides() -> None:
    server_guide = _read(SERVER_GUIDE)

    # ``cellucid serve`` prints the *absolute* path it resolved
    # (``anndata_server.py`` calls ``print_detail("Path", str(store.path))``),
    # so a photograph of this command necessarily bakes the capture host's path
    # into the image pixels, where no scan can find it. The transcripts therefore
    # ship as verbatim ``text`` blocks, never as screenshots.
    assert "```{figure} ../../../_static/screenshots/server/pancreas-cli-serve.png" not in (
        server_guide
    )
    assert "```{figure} ../../../_static/screenshots/server/serve-anndata-cli.png" not in (
        server_guide
    )
    assert "/?anndata=true" in server_guide
    assert "/?source=remote" not in server_guide
    assert "$ cellucid serve pbmc_demo.h5ad \\" in server_guide
    assert "      Format: h5ad file" in server_guide

    # The step counts are the part of a transcript that goes stale first,
    # because a new reported phase renumbers every line after it, so every
    # numbered header is pinned here.
    for anndata_step in (
        "[1/4] Detecting format...",
        "[2/4] Reading metadata...",
        "[3/4] Checking launch choices...",
        "[4/4] Starting server...",
    ):
        assert anndata_step in server_guide
    assert "[1/3]" not in server_guide
    assert "[1/5]" not in server_guide
    assert "[5/5]" not in server_guide

    # Step 2 reads metadata only, and step 3 prints every launch choice it
    # accepted, so the page quotes the lines a reader will be staring at.
    for step_line in (
        "      Cells: 4,096",
        "      obsm: X_umap",
        "      obsp: connectivities",
        "      embedding2d: {'key': 'X_umap', 'columns': [1, 0]}",
        "      knn: connectivities",
        "Consolidated metadata: built in memory (.zmetadata)",
        "Serving: a decodable copy, written once to a temporary .h5ad",
    ):
        assert step_line in server_guide

    # A bare ``X_umap`` fills the slot of its own column count, a named slot
    # wins, and only an object the default rule cannot read is guessed.
    # Describing that as anything other than a rule about width invites the
    # reader to rename a key that does not need renaming.
    normalized_server_guide = " ".join(server_guide.split())
    assert "fills the slot of its own column count" in normalized_server_guide
    assert "`--embedding-1d`, `--embedding-2d`, and `--embedding-3d`" in normalized_server_guide
    assert "deterministic best guess" in normalized_server_guide
    assert "Cellucid found nothing to plot." in server_guide
    assert "has shape (4096, 10)" not in server_guide

    # The viewer reads the store itself: the server publishes its bytes and
    # the launch choices, and none of the removed server-side work comes back.
    assert "/_cellucid/anndata/<dataset_id>.h5ad" in server_guide
    assert '"knn": "connectivities"' in server_guide
    assert "--knn connectivities" in server_guide
    assert 'knn_key="connectivities"' in server_guide
    for removed in (
        "serve_connectivity",
        "--connectivity",
        "--vector-fields",
        "serve_vector_fields",
        "backed",
        "_cellucid/mapping",
        "AnnDataAdapter",
    ):
        assert removed not in server_guide

    # A wildcard bind prints an origin a browser can open, never the wildcard.
    assert "Bound to every network interface. From another machine:" in server_guide
    assert "Local URL:    http://0.0.0.0" not in server_guide
    assert "Viewer URL:   http://0.0.0.0" not in server_guide
    assert "Machine URL:  http://compute-node-42:8765" in server_guide
    assert "ssh -N -L 8765:compute-node-42:8765 you@login-node" in server_guide
    assert (
        "{doc}`../../python_package/d_viewing_apis/"
        "12_remote_servers_ssh_tunneling_and_cloud`"
    ) in server_guide

    # ``cellucid serve`` no longer answers every ImportError with pip install.
    assert "no package installs a standard-library module" in normalized_server_guide
    assert "pip install --upgrade <name>" in server_guide

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
    picker = _read(DATA_LOADING / "03_browser_file_picker_tutorial.md")
    assert "**Open…** under **File or folder**, then **File…**" in overview
    assert "Browser **Folder** picker" not in overview
    assert "I chose Zarr folder…, but nothing happens" in picker


def test_python_and_r_sharing_navigation_points_to_the_reference_repository() -> None:
    paths = (
        DOCS_ROOT
        / "user_guide"
        / "python_package"
        / "d_viewing_apis"
        / "03_choose_your_workflow_decision_tree.md",
        DOCS_ROOT
        / "user_guide"
        / "r_package"
        / "d_viewing_loading"
        / "02_host_stores_for_sharing.md",
    )
    for path in paths:
        assert "11_custom_dataset_repository" in _read(path), path

    r_guide = _read(paths[-1])
    assert REFERENCE_VALUE in r_guide
    assert "`%||%`" not in r_guide
    assert "if (!file.exists(ident_path)) next" not in r_guide
