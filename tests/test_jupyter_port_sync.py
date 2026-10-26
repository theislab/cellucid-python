from __future__ import annotations

import html as html_module

import pytest


class _DummyServer:
    """Stands in for the AnnData server: records the bind, serves nothing."""

    port_offset = 0

    def __init__(self, **kwargs: object) -> None:
        self.host = kwargs["host"]
        self.port = int(kwargs["port"] or 8765) + self.port_offset

    @property
    def url(self) -> str:
        return f"http://{self.host}:{self.port}"

    def start_background(self) -> None:
        return None

    def is_running(self) -> bool:
        return True

    def stop(self) -> None:
        return None


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
        html = viewer._generate_viewer_html()
        assert html_module.escape(viewer.viewer_url, quote=True) in html
        assert 'src="about:blank"' not in html
        assert "fetch(" not in html
        assert "iframe.src =" not in html
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
        html = viewer._generate_viewer_html()
        assert explicit in html
        assert "ignored.example" not in html
        assert "127.0.0.1" not in html
        assert 'src="about:blank"' not in html
        assert "fetch(" not in html
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
