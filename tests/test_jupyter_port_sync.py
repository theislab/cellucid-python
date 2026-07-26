from __future__ import annotations

import html as html_module
from pathlib import Path

import pytest


def test_jupyter_viewer_port_sync(monkeypatch, tmp_path: Path) -> None:
    """
    If the bound server reports a different port, the viewer must use that exact port.
    """
    import cellucid.jupyter as jupyter

    class DummyServer:
        def __init__(
            self,
            data_dir: Path,
            port: int,
            host: str,
            open_browser: bool,
            quiet: bool,
            **_kwargs,
        ):
            self.data_dir = data_dir
            self.host = host
            self.port = int(port) + 1
            self.open_browser = open_browser
            self.quiet = quiet

        @property
        def url(self) -> str:
            return f"http://{self.host}:{self.port}"

        def start_background(self) -> None:
            return None

        def is_running(self) -> bool:
            return True

        def stop(self) -> None:
            return None

    monkeypatch.setattr(jupyter, "CellucidServer", DummyServer)

    viewer = jupyter.CellucidViewer(data_dir=tmp_path, port=8765, auto_open=False)
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


def test_jupyter_uses_only_the_explicit_browser_server_url(
    monkeypatch,
    tmp_path: Path,
) -> None:
    import cellucid.jupyter as jupyter

    class DummyServer:
        def __init__(self, **kwargs):
            self.host = kwargs["host"]
            self.port = kwargs["port"]

        @property
        def url(self) -> str:
            return f"http://{self.host}:{self.port}"

        def start_background(self) -> None:
            return None

        def is_running(self) -> bool:
            return True

        def stop(self) -> None:
            return None

    monkeypatch.setattr(jupyter, "CellucidServer", DummyServer)
    monkeypatch.setenv("CELLUCID_CLIENT_SERVER_URL", "https://ignored.example")

    explicit = "https://notebook.example/user/alice/gateway/8765"
    viewer = jupyter.CellucidViewer(
        data_dir=tmp_path,
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
    tmp_path: Path,
) -> None:
    import cellucid.jupyter as jupyter

    with pytest.raises((TypeError, ValueError), match="client_server_url"):
        jupyter.CellucidViewer(
            data_dir=tmp_path,
            auto_open=False,
            client_server_url=value,
        )
