"""The viewer token never reaches what is logged or pasted.

The address that carries it is said once, to the person, where they are told
to open it; an error, a debug report, a request log line and a ``repr`` say
the address without it (:func:`cellucid._server_base.public_viewer_url`, as
cellucid-r's ``.public_viewer_url()`` and ``test-viewer-token-in-output.R``).
"""

from __future__ import annotations

import ast
import http.client
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from cellucid._server_base import public_viewer_url
from cellucid.anndata_server import AnnDataServer

SOURCE = Path(__file__).resolve().parents[1] / "src" / "cellucid"

#: The places that say the tokened address to the person, to open it.
SAID_TO_OPEN = {
    # The banner of cellucid serve / serve_anndata(): "Viewer URL: ...".
    ("_server_base.py", "print_server_banner"),
    # AnnDataViewer.display() outside a notebook, and its notebook link.
    ("jupyter/_viewer.py", "display"),
    ("jupyter/_viewer.py", "_repr_html_"),
    # cellucid mcp --open: on stderr, the agent app's console.
    ("mcp.py", "_open_viewer"),
}


@pytest.mark.parametrize(
    ("text", "public"),
    [
        (
            "http://127.0.0.1:8765/?viewer=true&viewerId=ab&viewerToken=cd&anndata=true",
            "http://127.0.0.1:8765/?viewer=true&viewerId=ab&anndata=true",
        ),
        ("a?viewerToken=X&b=1", "a?b=1"),
        ("a?b=1&viewerToken=X", "a?b=1"),
        ("a?viewerToken=X", "a"),
        ("a?b=1", "a?b=1"),
        ("x?viewerToken=1#frag", "x#frag"),
        (
            "GET /?viewer=true&viewerToken=ff&anndata=true HTTP/1.1",
            "GET /?viewer=true&anndata=true HTTP/1.1",
        ),
        ("open a?viewerToken=1 or b?c=2&viewerToken=3", "open a or b?c=2"),
    ],
)
def test_the_public_address_drops_the_token_wherever_it_stands_as_r_does(text: str, public: str) -> None:
    assert public_viewer_url(text) == public


def _interpolated_viewer_urls(tree: ast.AST) -> list[tuple[str, int]]:
    """Each function (by name) that writes a ``viewer_url`` into text not
    through ``public_viewer_url``, with the line."""
    found: list[tuple[str, int]] = []

    def visit(node: ast.AST, function: str, stripped: bool) -> None:
        if isinstance(node, ast.FunctionDef | ast.AsyncFunctionDef):
            function = node.name
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "public_viewer_url":
            stripped = True
        if isinstance(node, ast.FormattedValue) and not stripped:
            for inner in ast.walk(node.value):
                named = inner.attr if isinstance(inner, ast.Attribute) else getattr(inner, "id", None)
                if named == "viewer_url" and not any(
                    isinstance(call, ast.Call) and getattr(call.func, "id", None) == "public_viewer_url"
                    for call in ast.walk(node.value)
                ):
                    found.append((function, node.lineno))
                    break
        for child in ast.iter_child_nodes(node):
            visit(child, function, stripped)

    visit(tree, "<module>", False)
    return found


def test_only_the_sentences_that_tell_the_person_what_to_open_say_the_tokened_address() -> None:
    unclassified = []
    for path in sorted(SOURCE.rglob("*.py")):
        relative = path.relative_to(SOURCE).as_posix()
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for function, line in _interpolated_viewer_urls(tree):
            if (relative, function) not in SAID_TO_OPEN:
                unclassified.append(f"{relative}:{line} ({function})")
    assert unclassified == [], "say public_viewer_url(...) there, or classify the sentence above"


def test_the_check_finds_an_address_written_into_an_error() -> None:
    tree = ast.parse(
        "def start(self):\n"
        "    raise RuntimeError(f'Could not open {self.viewer_url}')\n"
        "def fine(self):\n"
        "    raise RuntimeError(f'Could not open {public_viewer_url(self.viewer_url)}')\n"
    )
    assert _interpolated_viewer_urls(tree) == [("start", 2)]


def test_the_request_log_names_the_page_without_its_token(caplog: pytest.LogCaptureFixture) -> None:
    rng = np.random.default_rng(0)
    adata = AnnData(
        X=rng.random((4, 2)).astype(np.float32),
        obs=pd.DataFrame(index=[f"c{i}" for i in range(4)]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    adata.obsm["X_umap"] = rng.random((4, 2)).astype(np.float32)
    server = AnnDataServer(
        adata, host="127.0.0.1", port=0, quiet=True, serve_web_ui=False, dataset_id="t", dataset_name="T"
    )
    server.viewer_channel = True
    server.start_background()
    try:
        url = server.viewer_url
        token = url.split("viewerToken=", 1)[1].split("&", 1)[0]
        with caplog.at_level(logging.DEBUG, logger="cellucid.anndata_server"):
            connection = http.client.HTTPConnection("127.0.0.1", server.port, timeout=10)
            try:
                connection.request("GET", url.split(f":{server.port}", 1)[1])
                connection.getresponse().read()
            finally:
                connection.close()
        logged = "\n".join(record.getMessage() for record in caplog.records)
        assert "GET /?viewer=true" in logged
        assert token not in logged
    finally:
        server.stop()
