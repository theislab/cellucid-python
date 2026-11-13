"""What ``cellucid serve`` prints reaches the person at once, whatever its
standard output is.

Python block-buffers standard output when it is not a terminal: under
``cellucid serve data.h5ad --no-browser > serve.out`` (a process manager,
``nohup``, an IDE's run panel) a Viewer URL left in the buffer while the
server runs is one the person cannot open, so ``console_print`` flushes every
line it writes.
"""

from __future__ import annotations

import queue
import re
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

from _repository_source import CHILD_PROCESS_FENCE, repository_source_environment

_TOKEN = re.compile(r"viewerToken=[^&\s]+")


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", 0))
        return int(probe.getsockname()[1])


def _tiny_h5ad(path: Path) -> Path:
    import anndata as ad
    import numpy as np

    ad.AnnData(
        X=np.ones((3, 2), dtype=np.float32),
        obsm={"X_umap": np.zeros((3, 2), dtype=np.float32)},
    ).write_h5ad(path)
    return path


def test_cellucid_serve_says_its_viewer_url_into_a_pipe_while_it_runs(tmp_path: Path) -> None:
    environment = repository_source_environment(CELLUCID_AGENTS="off")
    environment.pop("PYTHONUNBUFFERED", None)
    command = [
        sys.executable,
        "-c",
        "import sys; from cellucid.cli import main; sys.exit(main(sys.argv[1:]))",
        "serve",
        str(_tiny_h5ad(tmp_path / "tiny.h5ad")),
        "--no-browser",
        "--no-web-ui",
        "--port",
        str(_free_port()),
    ]
    # Standard output is a pipe, as in `cellucid serve … > serve.out`.
    server = subprocess.Popen(command, env=environment, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    lines: queue.Queue[str | None] = queue.Queue()

    def read() -> None:
        assert server.stdout is not None
        for line in server.stdout:
            lines.put(line)
        lines.put(None)

    threading.Thread(target=read, daemon=True).start()
    said: list[str] = []
    viewer_url = False
    try:
        deadline = time.monotonic() + CHILD_PROCESS_FENCE
        while not viewer_url and time.monotonic() < deadline:
            try:
                line = lines.get(timeout=1)
            except queue.Empty:
                continue
            if line is None:
                break
            said.append(_TOKEN.sub("viewerToken=…", line.rstrip()))
            viewer_url = line.strip().startswith("Viewer URL:")
        # Read while the server still runs: nothing it said waited for it to stop.
        assert server.poll() is None, "\n".join(said)
    finally:
        server.terminate()
        try:
            server.wait(timeout=30)
        except subprocess.TimeoutExpired:
            server.kill()
            server.wait()
    # The line is never shown: the viewer's address carries its token.
    assert viewer_url, "no Viewer URL line before the server stopped:\n" + "\n".join(said)
    assert any(line.strip().startswith("Local URL:") for line in said)
