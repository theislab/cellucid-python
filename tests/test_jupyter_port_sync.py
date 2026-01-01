from __future__ import annotations

from pathlib import Path


def test_jupyter_viewer_port_sync(monkeypatch, tmp_path: Path) -> None:
    """
    If the underlying server bumps the port (e.g. requested port unavailable),
    the viewer should keep `viewer.port` in sync so notebook proxy URLs are correct.
    """
    import cellucid.jupyter as jupyter

    class DummyServer:
        def __init__(self, data_dir: Path, port: int, host: str, open_browser: bool, quiet: bool):
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
        assert "var port = 8766;" in html
    finally:
        viewer.stop()

