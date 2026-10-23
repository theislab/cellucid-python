"""A notebook viewer that never starts a data server.

Hooks, inbound events, outbound commands and the debug report are the
viewer's own work: testing them needs an :class:`AnnDataViewer` set up exactly
as its constructor sets one up, short of starting the
:class:`~cellucid.anndata_server.AnnDataServer` it would read the data through.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from cellucid._server_base import CELLUCID_WEB_URL
from cellucid.jupyter import AnnDataViewer


def serverless_viewer(
    *,
    port: int | None = None,
    height: int = 600,
    client_server_url: str | None = None,
    web_source_url: str = CELLUCID_WEB_URL,
    web_cache_dir: str | Path | None = None,
    allowed_hosts: Sequence[str] | None = None,
) -> AnnDataViewer:
    """Return an :class:`AnnDataViewer` with every part but its server."""
    viewer = object.__new__(AnnDataViewer)
    viewer._init_state(
        port=port,
        height=height,
        client_server_url=client_server_url,
        web_source_url=web_source_url,
        web_cache_dir=web_cache_dir,
        allowed_hosts=allowed_hosts,
    )
    return viewer
