"""The round-trip tutorials: seven Python notebooks and their seven R twins.

They are published on the docs site and open in Google Colab from their
badges, so each must be a clean notebook, link its twin, open in Colab at its
own path, download only what the repository holds, and never wait for the
viewer inside a cell (a notebook widget hears the viewer only between cells). Every Python call is also checked against the
current signatures by ``test_docs_current_api_contract.py``; the R calls are
checked against the R package's exports here.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from _web_root import sibling_file

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DOCS = REPOSITORY_ROOT / "docs" / "user_guide"
PYTHON_FOLDER = DOCS / "python_package" / "f_notebooks_tutorials"
R_FOLDER = DOCS / "r_package" / "f_notebooks_tutorials"
GALLERY = PYTHON_FOLDER / "06_cellucid_in_your_pipeline_python_and_r.md"
R_INDEX = R_FOLDER / "index.md"
R_REFERENCE = DOCS / "r_package" / "g_api_reference_coverage" / "01_public_functions_and_classes.md"
COLAB = "https://colab.research.google.com/github/theislab/cellucid-python/blob/main/"
RAW = "https://raw.githubusercontent.com/theislab/cellucid-python/main/"

#: Twin by twin: (Python notebook, R notebook).
TWINS = (
    ("41_show_a_dataset_in_any_notebook", "07_view_in_a_notebook"),
    ("42_pull_everything_back", "08_pull_everything_back"),
    ("43_push_clusters_refine_and_pull", "09_push_clusters_refine_and_pull"),
    ("44_round_trip_without_a_notebook", "10_round_trip_without_a_notebook"),
    ("45_gene_sets_round_trip", "11_gene_sets_round_trip"),
    ("46_export_an_h5ad_subset_and_reanalyze", "12_subset_and_reanalyze"),
    ("47_subsampled_views", "13_subsampled_views"),
)

#: Calls that wait for the viewer's answer: in a notebook widget they cannot
#: be answered inside the cell that waits.
PYTHON_WAITS = ("wait_for_ready(", "wait_for_event(", "get_session_bundle(")
R_WAITS = ("cellucid_wait_ready(", "cellucid_session(")


def _notebook(path: Path) -> dict:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    assert notebook["nbformat"] == 4, path
    ids = []
    for cell in notebook["cells"]:
        assert isinstance(cell.get("id"), str) and cell["id"], path
        ids.append(cell["id"])
        assert isinstance(cell["source"], list), path
        if cell["cell_type"] == "code":
            assert cell["execution_count"] is None and cell["outputs"] == [], (path, cell["id"])
    assert len(ids) == len(set(ids)), path
    return notebook


def _text(notebook: dict, kind: str | None = None) -> str:
    return "\n".join(
        "".join(cell["source"]) for cell in notebook["cells"] if kind is None or cell["cell_type"] == kind
    )


@pytest.mark.parametrize(("python", "r"), TWINS)
def test_each_tutorial_has_its_twin_its_colab_badge_and_no_waiting_cell(python: str, r: str) -> None:
    py_notebook = _notebook(PYTHON_FOLDER / f"{python}.ipynb")
    r_notebook = _notebook(R_FOLDER / f"{r}.ipynb")
    assert py_notebook["metadata"]["kernelspec"]["language"] == "python"
    assert r_notebook["metadata"]["kernelspec"] == {"display_name": "R", "language": "R", "name": "ir"}

    for notebook, folder, name, twin in (
        (py_notebook, PYTHON_FOLDER, python, f"r_package/f_notebooks_tutorials/{r}.ipynb"),
        (r_notebook, R_FOLDER, r, f"python_package/f_notebooks_tutorials/{python}.ipynb"),
    ):
        markdown = _text(notebook, "markdown")
        own = f"{COLAB}{folder.relative_to(REPOSITORY_ROOT).as_posix()}/{name}.ipynb"
        assert own in markdown, f"{name}: no Colab badge for its own path"
        assert twin in markdown, f"{name}: does not link its twin {twin}"
    py_code = _text(py_notebook, "code")
    r_code = _text(r_notebook, "code")
    for call in PYTHON_WAITS:
        assert call not in py_code, (python, call)
    for call in R_WAITS:
        assert call not in r_code, (r, call)

    gallery = GALLERY.read_text(encoding="utf-8")
    r_index = R_INDEX.read_text(encoding="utf-8")
    assert f"{{doc}}`{python}`" in gallery and f"r_package/f_notebooks_tutorials/{r}`" in gallery
    assert f"{{doc}}`{r}`" in r_index and f"python_package/f_notebooks_tutorials/{python}`" in r_index


def test_the_python_tutorials_download_only_what_the_repository_holds() -> None:
    for path in sorted(PYTHON_FOLDER.glob("4[0-9]_*.ipynb")):
        code = _text(_notebook(path), "code")
        for url in re.findall(r'"(https://raw\.githubusercontent\.com/[^"]+)"', code):
            assert url.startswith(RAW), (path.name, url)
        if RAW in code:
            folder = re.search(rf'"{re.escape(RAW)}"\s*"([^"]+)"', code)
            assert folder is not None, path.name
            name = re.search(r'NAME = "([^"]+)"', code)
            assert name is not None, path.name
            assert (REPOSITORY_ROOT / folder.group(1) / name.group(1)).is_file(), path.name


def _r_calls() -> list[tuple[str, str]]:
    """Each R tutorial and every cellucid function it calls."""
    return [
        (path.name, name)
        for path in sorted(R_FOLDER.glob("*.ipynb"))
        for name in sorted(set(re.findall(r"\b(cellucid_\w+)\(", _text(_notebook(path), "code"))))
    ]


def test_every_r_call_in_the_tutorials_is_documented() -> None:
    reference = R_REFERENCE.read_text(encoding="utf-8")
    for notebook, name in _r_calls():
        assert f"`{name}()`" in reference, f"{notebook}: {name}() is not in the R API reference"


def test_every_r_call_in_the_tutorials_is_exported_by_cellucid_r() -> None:
    namespace = sibling_file("cellucid-r", "NAMESPACE").read_text(encoding="utf-8")
    exported = set(re.findall(r"^export\((\w+)\)", namespace, re.MULTILINE))
    for notebook, name in _r_calls():
        assert name in exported, f"{notebook}: {name}() is not exported by cellucid-r"
