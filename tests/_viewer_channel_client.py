"""A fake viewer for the viewer channel: a plain-socket WebSocket client.

It speaks exactly what a browser does on ``/_cellucid/ws`` (masked client
frames, the RFC 6455 handshake with an ``Origin``) so the tests exercise the
server's real framing, fragmentation and close handling, not a library's.
"""

from __future__ import annotations

import base64
import contextlib
import json
import os
import socket
import struct
from typing import Any

from cellucid._websocket import accept_key

TEXT = 0x1
BINARY = 0x2
CLOSE = 0x8
PING = 0x9
PONG = 0xA


class HandshakeRefused(Exception):
    def __init__(self, status: int, body: bytes) -> None:
        super().__init__(f"HTTP {status}: {body!r}")
        self.status = status
        self.body = body


def binary_event(header: dict[str, Any], body: bytes) -> bytes:
    """One binary event frame payload: uint32 LE header length, JSON header, body."""
    encoded = json.dumps(header, separators=(",", ":")).encode("utf-8")
    return struct.pack("<I", len(encoded)) + encoded + body


class FakeViewer:
    def __init__(
        self,
        port: int,
        *,
        origin: str | None = None,
        host: str | None = None,
        path: str = "/_cellucid/ws",
        timeout: float = 5.0,
    ) -> None:
        self.sock = socket.create_connection(("127.0.0.1", port), timeout=timeout)
        key = base64.b64encode(os.urandom(16)).decode("ascii")
        host = host or f"127.0.0.1:{port}"
        origin = f"http://127.0.0.1:{port}" if origin is None else origin
        lines = [
            f"GET {path} HTTP/1.1",
            f"Host: {host}",
            "Upgrade: websocket",
            "Connection: Upgrade",
            f"Sec-WebSocket-Key: {key}",
            "Sec-WebSocket-Version: 13",
        ]
        if origin:
            lines.append(f"Origin: {origin}")
        self.sock.sendall(("\r\n".join(lines) + "\r\n\r\n").encode("ascii"))
        self._buffer = b""
        head = self._read_until(b"\r\n\r\n")
        status_line, *header_lines = head.decode("latin-1").split("\r\n")
        status = int(status_line.split()[1])
        headers = {
            name.strip().lower(): value.strip()
            for name, value in (line.split(":", 1) for line in header_lines if line)
        }
        if status != 101:
            length = int(headers.get("content-length", "0"))
            body = self._read_exact(length)
            self.sock.close()
            raise HandshakeRefused(status, body)
        assert headers["sec-websocket-accept"] == accept_key(key)

    # -- raw io --------------------------------------------------------------

    def _read_until(self, marker: bytes) -> bytes:
        while marker not in self._buffer:
            chunk = self.sock.recv(65536)
            if not chunk:
                raise ConnectionError("closed during handshake")
            self._buffer += chunk
        head, _, self._buffer = self._buffer.partition(marker)
        return head

    def _read_exact(self, size: int) -> bytes:
        while len(self._buffer) < size:
            chunk = self.sock.recv(max(65536, size - len(self._buffer)))
            if not chunk:
                raise ConnectionError("closed")
            self._buffer += chunk
        data, self._buffer = self._buffer[:size], self._buffer[size:]
        return data

    # -- frames ---------------------------------------------------------------

    @staticmethod
    def frame(opcode: int, payload: bytes, *, fin: bool = True, mask: bytes | None = None) -> bytes:
        mask = os.urandom(4) if mask is None else mask
        length = len(payload)
        first = (0x80 if fin else 0) | opcode
        if length < 126:
            header = struct.pack("!BB", first, 0x80 | length)
        elif length < 1 << 16:
            header = struct.pack("!BBH", first, 0x80 | 126, length)
        else:
            header = struct.pack("!BBQ", first, 0x80 | 127, length)
        tiled = (mask * (length // 4 + 1))[:length]
        masked = (int.from_bytes(payload, "little") ^ int.from_bytes(tiled, "little")).to_bytes(
            length, "little"
        )
        return header + mask + masked

    def send_frame(self, opcode: int, payload: bytes, *, fin: bool = True) -> None:
        self.sock.sendall(self.frame(opcode, payload, fin=fin))

    def send_json(self, value: dict[str, Any]) -> None:
        self.send_frame(TEXT, json.dumps(value).encode("utf-8"))

    def send_binary(self, payload: bytes) -> None:
        self.send_frame(BINARY, payload)

    def hello(self, viewer_id: str, viewer_token: str) -> dict[str, Any]:
        self.send_json({"type": "hello", "viewerId": viewer_id, "viewerToken": viewer_token})
        return self.receive_json()

    def receive(self) -> tuple[int, bytes]:
        """The next frame the server sends; pings are answered and skipped."""
        while True:
            first, second = self._read_exact(2)
            assert not second & 0x80, "server frames are never masked"
            length = second & 0x7F
            if length == 126:
                (length,) = struct.unpack("!H", self._read_exact(2))
            elif length == 127:
                (length,) = struct.unpack("!Q", self._read_exact(8))
            payload = self._read_exact(length)
            opcode = first & 0x0F
            if opcode == PING:
                self.send_frame(PONG, payload)
                continue
            return opcode, payload

    def receive_json(self) -> dict[str, Any]:
        opcode, payload = self.receive()
        if opcode == CLOSE:
            raise ConnectionError(f"closed {struct.unpack('!H', payload[:2])[0]}")
        assert opcode == TEXT
        return json.loads(payload.decode("utf-8"))

    def receive_command(self) -> dict[str, Any]:
        """The next command, text or binary: a binary one as its header with
        its ``cells`` (an empty list for ``setVisibility`` of every cell), the
        header's length checked to put the cells on a 4-byte boundary."""
        opcode, payload = self.receive()
        if opcode == TEXT:
            return json.loads(payload.decode("utf-8"))
        assert opcode == BINARY, opcode
        (header_length,) = struct.unpack("<I", payload[:4])
        assert header_length % 4 == 0
        header = json.loads(payload[4 : 4 + header_length].decode("utf-8"))
        body = payload[4 + header_length :]
        cells = list(struct.unpack(f"<{len(body) // 4}I", body))
        return {**header, "cells": cells}

    def receive_close(self) -> tuple[int, str]:
        """The server's close frame, its code and reason.

        The close must be orderly (RFC 6455 section 7.1.1): the server shuts
        its write side down right after the frame, so the stream ends there.
        A reset instead -- the server closing its socket with input unread --
        is what Windows reports in place of a close frame the client has not
        read yet, so it fails here on every platform.
        """
        while True:
            opcode, payload = self.receive()
            if opcode == CLOSE:
                assert self._buffer == b"", "the server sent bytes after its close frame"
                assert self.sock.recv(1) == b"", "the server sent bytes after its close frame"
                return struct.unpack("!H", payload[:2])[0], payload[2:].decode("utf-8")

    def answer_close(self) -> None:
        """Answer the server's close as a browser does: a close frame, then the
        client's own end of the connection."""
        self.send_frame(CLOSE, struct.pack("!H", 1000))
        self.sock.shutdown(socket.SHUT_WR)

    def close(self) -> None:
        with contextlib.suppress(OSError):
            self.send_frame(CLOSE, struct.pack("!H", 1000))
        self.sock.close()
