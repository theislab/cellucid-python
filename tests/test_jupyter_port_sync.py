from __future__ import annotations

import html as html_module

import pytest
from _stand_in_server import StandInServer as _DummyServer


def test_jupyter_viewer_port_sync(monkeypatch) -> None:
    """If the bound server reports a different port, the viewer must use that exact port."""
    import cellucid.anndata_server as anndata_server
    import cellucid.jupyter as jupyter

    class BoundElsewhere(_DummyServer):
        port_offset = 1

    # ``AnnDataViewer._start_server`` imports ``AnnDataServer`` from its module
    # when it runs, so the patch names that module.
    monkeypatch.setattr(anndata_server, "AnnDataServer", BoundElsewhere)

    viewer = jupyter.AnnDataViewer("data.h5ad", port=8765, auto_open=False)
    try:
        assert viewer.port == 8766
        assert "8766" in viewer.viewer_url
        html = viewer._repr_html_()
        assert html_module.escape(viewer.viewer_url, quote=True) in html
        assert "<iframe" not in html
        assert "<script" not in html
    finally:
        viewer.stop()


def test_jupyter_uses_only_the_explicit_browser_server_url(monkeypatch) -> None:
    import cellucid.anndata_server as anndata_server
    import cellucid.jupyter as jupyter

    monkeypatch.setattr(anndata_server, "AnnDataServer", _DummyServer)
    monkeypatch.setenv("CELLUCID_CLIENT_SERVER_URL", "https://ignored.example")

    explicit = "https://notebook.example/user/alice/gateway/8765"
    viewer = jupyter.AnnDataViewer(
        "data.h5ad",
        port=8765,
        auto_open=False,
        client_server_url=explicit,
    )
    try:
        assert viewer.viewer_url.startswith(f"{explicit}/?")
        html = viewer._repr_html_()
        assert explicit in html
        assert "ignored.example" not in html
        assert "127.0.0.1" not in html
    finally:
        viewer.stop()


@pytest.mark.parametrize(
    "value",
    [
        "",
        "not-a-url",
        "ftp://example.com",
        "https://example.com/",
        "https://example.com:invalid",
        " https://example.com",
        "https://example.com/path?query=1",
        "https://example.com/path#fragment",
    ],
)
def test_jupyter_rejects_invalid_explicit_browser_server_urls(
    value: str,
    monkeypatch,
) -> None:
    import cellucid.anndata_server as anndata_server
    import cellucid.jupyter as jupyter

    monkeypatch.setattr(anndata_server, "AnnDataServer", _DummyServer)
    with pytest.raises((TypeError, ValueError), match="client_server_url"):
        jupyter.AnnDataViewer(
            "data.h5ad",
            auto_open=False,
            client_server_url=value,
        )


def test_the_stand_in_server_has_what_the_viewers_read() -> None:
    """Every attribute the viewers read from their server (``self._server.<name>``
    in ``jupyter/_viewer.py`` and ``_widget.py``) is on the one stand-in the
    tests use, so a new one is added there once, not missed in a copy."""
    import re
    from pathlib import Path

    import cellucid.jupyter as jupyter

    package = Path(jupyter.__file__).parent
    read = {
        name
        for module in ("_viewer.py", "_widget.py")
        for name in re.findall(r"self\._server\.([A-Za-z_]+)", (package / module).read_text(encoding="utf-8"))
    }
    assert read, "the viewers read nothing from their server: the pattern is stale"
    stand_in = _DummyServer(host="127.0.0.1", port=None)
    assert sorted(name for name in read if not hasattr(stand_in, name)) == []
