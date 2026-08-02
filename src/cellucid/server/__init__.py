"""
Cellucid Data Server

A lightweight HTTP server for serving pre-exported cellucid datasets.
Supports both local and remote access patterns:

1. Local mode: Run on your machine, open browser locally
2. SSH tunnel mode: Run on remote server, access via SSH port forwarding
3. Jupyter mode: Run alongside Jupyter, embed in notebooks

Usage:
    from cellucid import serve
    serve("/path/to/data", port=8765)

Or via CLI:
    cellucid serve /path/to/data --port 8765

The server provides:
- Static file serving for dataset files
- CORS headers for cross-origin access (needed for web viewer)
- Health check endpoint for connection validation
- A wire capability endpoint (``/_cellucid/protocol``) the viewer reads before
  it emits an event this ``cellucid`` may not accept

For serving AnnData directly (without pre-export), use:
    from cellucid import serve_anndata
    serve_anndata(
        "/path/to/data.h5ad",
        dataset_name="Example",
        dataset_id="example",
    )

The module is a package. Its public surface is the two names below plus the
request handler; everything else is split by responsibility across private
submodules:

``_state``
    The published-state sidecars an export directory may carry.
``_datasets``
    Discovery of the exported datasets a served root declares.
``_artifacts``
    The exact inventory of files one dataset is allowed to serve.
``_handler``
    :class:`CORSRequestHandler` and the ``cellucid.server`` logger.
``_server``
    :class:`CellucidServer` and :func:`serve`.
"""

from ._handler import CORSRequestHandler
from ._server import CellucidServer, serve

__all__ = ["CORSRequestHandler", "CellucidServer", "serve"]
