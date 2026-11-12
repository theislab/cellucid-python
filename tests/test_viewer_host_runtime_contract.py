"""What opened a viewer, as each entry point names it.

The viewer names its host in every sentence it says from the ``runtime`` of
``/_cellucid/info`` (the shared protocol's ``runtimes``): ending ``cellucid
mcp --open`` once froze the page with "The notebook or R session that opened
this viewer stopped". Each way this package opens a viewer names itself:
a notebook kernel ``jupyter``, a script or session ``python``, ``cellucid
serve`` ``cellucid-serve``, ``cellucid mcp --open`` ``cellucid-mcp``; and
the agents record it publishes carries the same name.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from _stand_in_server import StandInServer

PROTOCOL = json.loads((Path(__file__).parent / "fixtures" / "viewer_protocol.json").read_text(encoding="utf-8"))
PYTHON_RUNTIMES = {"jupyter", "python", "cellucid-serve", "cellucid-mcp"}


class _Server(StandInServer):
    """The stand-in server, recording what opened it."""

    made: list[_Server] = []

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.started_as: str | None = None
        _Server.made.append(self)

    def start(self) -> None:
        self.started_as = self.runtime

    def start_background(self) -> None:
        self.started_as = self.runtime


@pytest.fixture
def servers(monkeypatch: pytest.MonkeyPatch) -> list[_Server]:
    import cellucid.anndata_server as anndata_server

    _Server.made = []
    monkeypatch.setattr(anndata_server, "AnnDataServer", _Server)
    return _Server.made


@pytest.fixture
def published(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    import cellucid.jupyter._viewer as viewer_module

    runtimes: list[str] = []
    monkeypatch.setattr(viewer_module, "publish_viewer", lambda _viewer, _port, *, runtime: runtimes.append(runtime))
    return runtimes


def _in_notebook(monkeypatch: pytest.MonkeyPatch, in_jupyter: bool) -> None:
    import cellucid.jupyter._viewer as viewer_module

    detect = viewer_module._detect_jupyter_context

    def detected() -> dict[str, Any]:
        return {**detect(), "in_jupyter": in_jupyter}

    monkeypatch.setattr(viewer_module, "_detect_jupyter_context", detected)


def test_the_protocol_declares_every_runtime_this_package_names_as_python() -> None:
    runtimes = PROTOCOL["runtimes"]
    assert set(runtimes) == PYTHON_RUNTIMES | {"r", "r-mcp"}
    for runtime in PYTHON_RUNTIMES:
        assert runtimes[runtime]["language"] == "python", runtime


@pytest.mark.parametrize(("in_jupyter", "runtime"), [(False, "python"), (True, "jupyter")])
def test_an_anndata_viewer_names_the_python_that_opened_it(
    monkeypatch: pytest.MonkeyPatch, servers: list[_Server], published: list[str], in_jupyter: bool, runtime: str
) -> None:
    from cellucid import jupyter

    _in_notebook(monkeypatch, in_jupyter)
    viewer = jupyter.AnnDataViewer("data.h5ad", auto_open=False)
    try:
        assert servers[0].started_as == runtime
        assert published == [runtime]
        assert viewer._agent_description()["runtime"] == runtime
    finally:
        viewer.stop()


@pytest.mark.parametrize(("in_jupyter", "runtime"), [(False, "python"), (True, "jupyter")])
def test_serve_anndata_names_the_python_that_called_it(
    monkeypatch: pytest.MonkeyPatch, servers: list[_Server], in_jupyter: bool, runtime: str
) -> None:
    from cellucid import serve_anndata

    _in_notebook(monkeypatch, in_jupyter)
    serve_anndata("data.h5ad", open_browser=False, quiet=True)
    assert servers[0].viewer_channel is True
    assert servers[0].started_as == runtime


def test_cellucid_serve_names_itself(tmp_path: Path, servers: list[_Server]) -> None:
    import anndata as ad
    import numpy as np

    from cellucid.cli import main

    path = tmp_path / "tiny.h5ad"
    ad.AnnData(X=np.ones((3, 2), dtype=np.float32), obsm={"X_umap": np.zeros((3, 2), dtype=np.float32)}).write_h5ad(path)
    assert main(["serve", str(path), "--no-browser", "--quiet"]) == 0
    assert servers[0].viewer_channel is True
    assert servers[0].started_as == "cellucid-serve"


def test_cellucid_mcp_open_names_itself(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    servers: list[_Server],
    published: list[str],
) -> None:
    from cellucid import mcp

    _in_notebook(monkeypatch, True)  # even from a notebook's Python
    viewer = mcp._open_viewer("data.h5ad", port=None, open_browser=False, web_source_url=None, web_cache_dir=None)
    try:
        assert servers[0].started_as == "cellucid-mcp"
        assert published == ["cellucid-mcp"]
        assert viewer._agent_description()["runtime"] == "cellucid-mcp"
    finally:
        viewer.stop()
    # Its address went to stderr (read here, never shown: it carries the token).
    assert capsys.readouterr().err.startswith("Cellucid viewer: http://127.0.0.1:8765/")
