"""
The HTTP request handler for a prepared-data server.

:class:`CORSRequestHandler` answers only from the artifact inventory the
dataset declares, resolves every request path against the served root before
opening anything, and serves byte ranges for the payloads the viewer streams.
The shared CORS, web-asset proxy, and event/session-upload behaviour comes from
:class:`cellucid._server_base.CORSMixin`.

The module also owns the ``cellucid.server`` logger, which the server lifecycle
in :mod:`._server` shares.
"""

from __future__ import annotations

import logging
import os
import re
import stat
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler
from pathlib import Path
from urllib.parse import urlparse

from .._server_base import WEB_ASSET_INVENTORY_FILENAME, CORSMixin
from ._artifacts import _PreparedArtifact
from ._datasets import _list_exported_datasets, _read_exported_dataset_entry

logger = logging.getLogger("cellucid.server")


def _protocol_capability_document() -> dict[str, list[str]]:
    """Return the viewer wire capabilities this installation declares.

    The answer is read from :mod:`cellucid.jupyter._wire`, which is the module
    the validator itself reads, so the route cannot publish a capability the
    notebook would then reject. The import is deferred because
    ``cellucid/jupyter/_exported.py`` imports :class:`cellucid.server`
    at module level: naming ``_wire`` at the top of this module makes
    ``import cellucid.server`` circular.
    """
    from ..jupyter import _wire

    return _wire._protocol_capability_document()


def _parse_byte_range(value: str | None, size: int) -> tuple[int, int] | None:
    if value is None:
        return None
    if type(size) is not int or size < 0:
        raise ValueError("Artifact size must be a non-negative integer.")
    match = re.fullmatch(r"bytes=(\d{0,20})-(\d{0,20})", value)
    if match is None or (not match.group(1) and not match.group(2)):
        raise ValueError("Range must contain one exact byte interval.")
    if size == 0:
        raise ValueError("An empty artifact has no satisfiable byte range.")
    if match.group(1):
        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) else size - 1
        if start >= size or end < start:
            raise ValueError("Range is outside the artifact.")
        return start, min(end, size - 1)
    suffix_length = int(match.group(2))
    if suffix_length <= 0:
        raise ValueError("Range suffix length must be positive.")
    return max(0, size - suffix_length), size - 1


class CORSRequestHandler(CORSMixin, SimpleHTTPRequestHandler):
    """Serve only endpoints and immutable prepared artifacts in the active contract."""

    allow_caching = True  # Static files can be cached

    def __init__(
        self,
        *args,
        data_dir: Path,
        server_info: dict,
        datasets: list[dict[str, str]],
        artifact_inventory: dict[str, _PreparedArtifact],
        serve_web_ui: bool,
        web_cache_dir: Path,
        allowed_hosts: frozenset[str],
        **kwargs,
    ):
        self.data_dir = data_dir
        self.server_info = server_info
        self.datasets = datasets
        self.artifact_inventory = artifact_inventory
        self.serve_web_ui = serve_web_ui
        self.web_cache_dir = web_cache_dir
        # Read by CORSMixin.parse_request, which runs inside super().__init__().
        self.allowed_hosts = allowed_hosts
        # Must call super().__init__ last because it calls do_GET immediately
        super().__init__(*args, directory=str(data_dir), **kwargs)

    def end_headers(self):
        """Add CORS headers to every response."""
        self.add_cors_headers()
        super().end_headers()

    def do_POST(self):
        """Handle POST requests (events from frontend)."""
        if self.handle_event_post():
            return
        if self.handle_session_bundle_post():
            return
        # No other POST endpoints - return 404
        self.send_error_response(404, f"POST not supported for path: {self.path}")

    def do_GET(self):
        """Serve one exact GET endpoint or prepared artifact."""
        self._handle_read_request(head_only=False)

    def do_HEAD(self):
        """Serve GET metadata without a response body."""
        self._handle_read_request(head_only=True)

    def _handle_read_request(self, *, head_only: bool) -> None:
        path = self._canonical_request_path()
        if path is None:
            self.send_error_response(404, "Request path is not in the active contract")
            return

        if path == "/_cellucid/health":
            self.send_json(
                {
                    "status": "ok",
                    "type": "exported",
                    "version": self.server_info["version"],
                },
                head_only=head_only,
            )
            return

        if path == "/_cellucid/info":
            self.send_json(self.server_info, head_only=head_only)
            return

        if path == "/_cellucid/datasets":
            self.send_json({"datasets": self.datasets}, head_only=head_only)
            return

        if path == "/_cellucid/protocol":
            self.send_json(_protocol_capability_document(), head_only=head_only)
            return

        if self.serve_web_ui and self.handle_web_asset_get(path, head_only=head_only):
            return
        if self.serve_web_ui and (
            path == f"/{WEB_ASSET_INVENTORY_FILENAME}"
            or path == "/assets"
            or path.startswith("/assets/")
        ):
            self.send_error_response(404, "Web asset is not declared by the active build")
            return

        if path in {"/", "/index.html"}:
            self.send_error_response(503, "Cellucid viewer UI unavailable")
            return

        self._serve_prepared_artifact(path[1:], head_only=head_only)

    def _canonical_request_path(self) -> str | None:
        """Return one exact ASCII origin-form path without decoding aliases."""
        parsed = urlparse(self.path)
        exact_target = parsed.path
        if parsed.query:
            exact_target += f"?{parsed.query}"
        if (
            parsed.scheme
            or parsed.netloc
            or parsed.params
            or parsed.fragment
            or exact_target != self.path
            or (parsed.query and parsed.path != "/")
            or not parsed.path.startswith("/")
            or "\\" in parsed.path
            or "%" in parsed.path
        ):
            return None
        try:
            parsed.path.encode("ascii")
        except UnicodeEncodeError:
            return None
        if any(
            ord(character) < 33 or ord(character) == 127
            for character in parsed.path
            if character != "/"
        ):
            return None
        if parsed.path != "/":
            parts = parsed.path[1:].split("/")
            if any(part in {"", ".", ".."} for part in parts):
                return None
        return parsed.path

    @staticmethod
    def _artifact_metadata_matches(
        metadata: os.stat_result,
        artifact: _PreparedArtifact,
    ) -> bool:
        return (
            stat.S_ISREG(metadata.st_mode)
            and metadata.st_size == artifact.size
            and metadata.st_mtime_ns == artifact.mtime_ns
            and metadata.st_dev == artifact.device
            and metadata.st_ino == artifact.inode
        )

    def _open_prepared_artifact(self, artifact: _PreparedArtifact) -> int:
        """Open one unchanged regular artifact without following its final symlink."""
        try:
            if artifact.path.is_symlink() or artifact.path.resolve(strict=True) != artifact.path:
                raise OSError("Prepared artifact path changed.")
            before = artifact.path.stat()
            if not self._artifact_metadata_matches(before, artifact):
                raise OSError("Prepared artifact metadata changed.")

            flags = os.O_RDONLY | getattr(os, "O_BINARY", 0)
            nofollow = getattr(os, "O_NOFOLLOW", 0)
            if nofollow:
                flags |= nofollow
            descriptor = os.open(artifact.path, flags)
        except (OSError, RuntimeError):
            raise

        try:
            opened = os.fstat(descriptor)
            after = artifact.path.stat()
            if not self._artifact_metadata_matches(
                opened, artifact
            ) or not self._artifact_metadata_matches(after, artifact):
                raise OSError("Prepared artifact changed while it was opened.")
        except BaseException:
            os.close(descriptor)
            raise
        return descriptor

    def _serve_prepared_artifact(self, request_path: str, *, head_only: bool) -> None:
        artifact = self.artifact_inventory.get(request_path)
        if artifact is None:
            self.send_error_response(404, "Prepared artifact is not declared")
            return

        range_values = self._request_headers().get_all("Range", [])
        if len(range_values) > 1:
            self._send_range_error(artifact.size)
            return
        try:
            interval = _parse_byte_range(
                range_values[0] if range_values else None,
                artifact.size,
            )
        except ValueError:
            self._send_range_error(artifact.size)
            return

        try:
            descriptor = self._open_prepared_artifact(artifact)
        except OSError:
            self.send_error_response(
                HTTPStatus.CONFLICT,
                "Prepared artifact changed after server validation",
            )
            return

        start, end = interval if interval is not None else (0, artifact.size - 1)
        content_length = 0 if artifact.size == 0 else end - start + 1
        self.send_response(HTTPStatus.PARTIAL_CONTENT if interval is not None else HTTPStatus.OK)
        self.send_header("Content-Type", artifact.content_type)
        self.send_header("Content-Length", str(content_length))
        self.send_header("Accept-Ranges", "bytes")
        if interval is not None:
            self.send_header(
                "Content-Range",
                f"bytes {start}-{end}/{artifact.size}",
            )
        self.end_headers()

        try:
            if not head_only and content_length:
                os.lseek(descriptor, start, os.SEEK_SET)
                remaining = content_length
                while remaining:
                    chunk = os.read(descriptor, min(1024 * 1024, remaining))
                    if not chunk:
                        raise OSError("Prepared artifact ended during response.")
                    self._response_writer().write(chunk)
                    remaining -= len(chunk)
        finally:
            os.close(descriptor)

    def _send_range_error(self, artifact_size: int) -> None:
        body = b"Requested byte range is not satisfiable"
        self.send_response(HTTPStatus.REQUESTED_RANGE_NOT_SATISFIABLE)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Content-Range", f"bytes */{artifact_size}")
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        if self.command != "HEAD":
            self._response_writer().write(body)

    def _list_datasets(self) -> list[dict]:
        """Validate and return the prepared datasets under ``data_dir``."""
        return _list_exported_datasets(self.data_dir)

    def _is_dataset_dir(self, path: Path) -> bool:
        """Check if a directory is a valid cellucid dataset."""
        try:
            _read_exported_dataset_entry(path, public_path="/")
        except (OSError, TypeError, ValueError):
            return False
        return True

    def _get_dataset_identity_fields(self, path: Path) -> tuple[str, str]:
        """Return (dataset_id, dataset_name) for a dataset directory."""
        entry = _read_exported_dataset_entry(path, public_path="/")
        return entry["id"], entry["name"]

    def log_message(self, format: str, *args):
        """Override to use Python logging instead of stderr."""
        logger.debug("%s - %s", self.address_string(), format % args)
