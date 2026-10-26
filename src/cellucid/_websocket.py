"""
The server side of one WebSocket (RFC 6455) over an upgraded HTTP connection.

This module knows framing only: the opening handshake's accept key, reading
client frames (masked, possibly fragmented, with control frames interleaved),
and writing server frames. What the frames mean is
:mod:`cellucid._viewer_channel`'s business.

Reading is built for large binary messages: a message is read as a stream
(:class:`BinaryMessage`), its payload unmasked in place in the caller's buffer
(NumPy XOR over 32-bit words), so a selection of millions of cells or a session
bundle of hundreds of megabytes is never copied or held twice. Writing is
serialized by one lock, so any thread may send while the connection's own
thread reads.

Closing is the orderly close of RFC 6455 section 7.1.1: the close frame is
sent, the write side is shut down, and the connection's own thread reads and
discards what the peer still sends until it closes its side, for at most
:data:`CLOSE_DRAIN_SECONDS`, before the socket is closed (:meth:`WebSocket.finish`).
A socket closed with received bytes still unread is reset rather than shut
down, and on Windows the reset discards the close frame the peer has not read
yet, so the page would see an abnormal closure (1006) instead of the code and
reason it was sent.
"""

from __future__ import annotations

import base64
import contextlib
import hashlib
import math
import select
import socket
import struct
import sys
import threading
import time
from collections.abc import Callable
from typing import BinaryIO

import numpy as np

from ._plural import plural

_ACCEPT_GUID = b"258EAFA5-E914-47DA-95CA-C5AB0DC85B11"

OPCODE_CONTINUATION = 0x0
OPCODE_TEXT = 0x1
OPCODE_BINARY = 0x2
OPCODE_CLOSE = 0x8
OPCODE_PING = 0x9
OPCODE_PONG = 0xA

CLOSE_NORMAL = 1000
CLOSE_GOING_AWAY = 1001
CLOSE_PROTOCOL_ERROR = 1002
CLOSE_INVALID_DATA = 1007
CLOSE_POLICY_VIOLATION = 1008
CLOSE_MESSAGE_TOO_BIG = 1009

# Largest payload written in one piece with its header; a larger one is
# written after it, without being copied.
_COALESCE_BYTES = 64 * 1024

#: Longest the connection waits, after its close frame went out, for the peer
#: to close its side. A peer that is still sending (the rest of a frame the
#: server refused from its header) or that never answers is cut off then.
CLOSE_DRAIN_SECONDS = 5.0

#: How often a thread waiting for the peer's bytes checks whether another
#: thread closed the connection meanwhile.
_CLOSED_CHECK_SECONDS = 0.25

# Read size of the read-ahead and of the drain after a close.
_READ_AHEAD_BYTES = 64 * 1024

# Longest close reason a close frame can carry: a control frame's payload is
# at most 125 bytes, two of which are the code.
_CLOSE_REASON_BYTES = 123


class WebSocketClosed(ConnectionError):
    """The connection is closed: by the peer, by this side, or by the network."""

    def __init__(self, code: int, reason: str = "") -> None:
        super().__init__(f"WebSocket closed ({code}){': ' + reason if reason else ''}")
        self.code = code
        self.reason = reason


class WebSocketProtocolError(ValueError):
    """The peer broke the framing rules; the connection is closed with ``code``."""

    def __init__(self, code: int, reason: str) -> None:
        super().__init__(reason)
        self.code = code
        self.reason = reason


def accept_key(client_key: str) -> str:
    """The ``Sec-WebSocket-Accept`` answer to one ``Sec-WebSocket-Key``."""
    digest = hashlib.sha1(client_key.encode("ascii") + _ACCEPT_GUID).digest()
    return base64.b64encode(digest).decode("ascii")


def require_client_key(value: object) -> str:
    """Return a valid ``Sec-WebSocket-Key`` (16 random bytes, base64), else raise."""
    if type(value) is not str or len(value) != 24:
        raise ValueError("Sec-WebSocket-Key must be 16 bytes in base64")
    try:
        decoded = base64.b64decode(value, validate=True)
    except ValueError as error:
        raise ValueError("Sec-WebSocket-Key must be 16 bytes in base64") from error
    if len(decoded) != 16:
        raise ValueError("Sec-WebSocket-Key must be 16 bytes in base64")
    return value


def unmask_in_place(buffer: memoryview, mask: bytes, offset: int) -> None:
    """XOR ``buffer`` with the frame ``mask``, ``offset`` bytes into the payload.

    The bulk runs as one NumPy XOR over 32-bit words of the buffer itself; the
    at most three bytes before the first word boundary, and the tail, are done
    byte by byte. Nothing is allocated in proportion to the buffer.
    """
    length = len(buffer)
    if length == 0:
        return
    data = np.frombuffer(buffer, dtype=np.uint8)
    # Bytes up to the first payload offset that is a multiple of four, so the
    # word-wide XOR below always starts on mask byte 0.
    head = (-offset) % 4
    head = min(head, length)
    for index in range(head):
        data[index] ^= mask[(offset + index) % 4]
    body = length - head
    words = body // 4
    if words:
        start = head
        # A view of the buffer as uint32 needs 4-byte alignment of its start
        # address; a misaligned buffer is XORed through a byte-tiled mask view.
        segment = data[start : start + words * 4]
        if segment.ctypes.data % 4 == 0:
            np.bitwise_xor(
                segment.view(np.uint32),
                np.frombuffer(mask, dtype=np.uint32)[0],
                out=segment.view(np.uint32),
            )
        else:
            tiled = np.frombuffer(mask * 1024, dtype=np.uint8)
            position = 0
            while position < segment.size:
                count = min(tiled.size, segment.size - position)
                np.bitwise_xor(
                    segment[position : position + count],
                    tiled[:count],
                    out=segment[position : position + count],
                )
                position += count
    tail_start = head + words * 4
    for index in range(tail_start, length):
        data[index] ^= mask[(offset + index) % 4]


def _readiness_poller(connection: socket.socket) -> Callable[[float], bool]:
    """A wait for ``connection`` to have bytes or an end to read, in seconds;
    whether it does. ``poll`` where the platform has it (no descriptor limit),
    else ``select`` (Windows)."""
    if sys.platform != "win32":
        poller = select.poll()
        poller.register(connection, select.POLLIN)
        return lambda seconds: bool(poller.poll(max(0, math.ceil(seconds * 1000))))
    return lambda seconds: bool(select.select([connection], [], [], max(0.0, seconds))[0])


class _Frame:
    __slots__ = ("fin", "mask", "opcode", "remaining")

    def __init__(self, fin: bool, opcode: int, remaining: int, mask: bytes) -> None:
        self.fin = fin
        self.opcode = opcode
        self.remaining = remaining
        self.mask = mask


class BinaryMessage:
    """One binary message, read as a stream across its frames.

    :meth:`readinto` fills a caller's buffer with unmasked payload bytes and
    returns how many it wrote (0 at the end of the message). Control frames
    between fragments are answered transparently. A message that is not read
    to its end must be :meth:`drain`-ed before the next one is received.
    """

    __slots__ = ("_frame", "_offset", "_socket", "done")

    def __init__(self, socket_: WebSocket, frame: _Frame) -> None:
        self._socket = socket_
        self._frame = frame
        self._offset = 0
        self.done = frame.remaining == 0 and frame.fin

    def readinto(self, buffer: memoryview) -> int:
        """Fill ``buffer`` from the message; return the bytes written, 0 at its end."""
        written = 0
        while written < len(buffer) and not self.done:
            if self._frame.remaining == 0:
                self._frame = self._socket._next_data_frame(continuation=True)
                self._offset = 0
                if self._frame.remaining == 0:
                    self.done = self._frame.fin
                continue
            count = min(len(buffer) - written, self._frame.remaining)
            target = buffer[written : written + count]
            self._socket._read_into(target)
            unmask_in_place(target, self._frame.mask, self._offset)
            self._offset += count
            self._frame.remaining -= count
            written += count
            if self._frame.remaining == 0 and self._frame.fin:
                self.done = True
        return written

    def read_exact(self, size: int) -> bytearray:
        """Read exactly ``size`` bytes of the message, or raise if it ends first."""
        result = bytearray(size)
        filled = 0
        with memoryview(result) as view:
            while filled < size:
                count = self.readinto(view[filled:])
                if count == 0:
                    raise WebSocketProtocolError(
                        CLOSE_POLICY_VIOLATION, "Binary message ended early"
                    )
                filled += count
        return result

    def read_rest(self, initial: int = 0) -> bytearray:
        """Read the rest of the message into one buffer, sized from its frames."""
        result = bytearray(max(initial, self._frame.remaining))
        filled = 0
        while True:
            if filled == len(result):
                if self.done:
                    break
                result.extend(bytes(max(self._frame.remaining, 64 * 1024)))
            with memoryview(result) as whole, whole[filled:] as view:
                count = self.readinto(view)
            if count == 0:
                break
            filled += count
        del result[filled:]
        return result

    def drain(self) -> None:
        """Read and discard what is left of the message."""
        scratch = memoryview(bytearray(256 * 1024))
        while self.readinto(scratch):
            pass


class WebSocket:
    """The server end of one upgraded connection.

    One thread, the connection's own, reads it and calls :meth:`finish` when
    it is done with it; any thread may write and close.
    """

    def __init__(
        self,
        connection: socket.socket,
        writer: BinaryIO,
        *,
        received: bytes,
        max_text_bytes: int,
    ) -> None:
        """``received`` is what the peer sent that was read along with its
        handshake, before the connection was handed over."""
        self._connection = connection
        self._writer = writer
        self._max_text_bytes = max_text_bytes
        # Bytes read from the socket and not yet consumed.
        self._received = bytearray(received)
        self._scratch = memoryview(bytearray(_READ_AHEAD_BYTES))
        self._poller = _readiness_poller(connection)
        self._write_lock = threading.Lock()
        self._closed = threading.Event()
        # When the close frame went out: the end of the drain that follows it.
        self._drain_deadline: float | None = None
        self.close_code: int | None = None

    @property
    def closed(self) -> bool:
        return self._closed.is_set()

    # ------------------------------------------------------------------
    # Writing
    # ------------------------------------------------------------------

    @staticmethod
    def _frame_header(opcode: int, length: int) -> bytes:
        if length < 126:
            return struct.pack("!BB", 0x80 | opcode, length)
        if length < 1 << 16:
            return struct.pack("!BBH", 0x80 | opcode, 126, length)
        return struct.pack("!BBQ", 0x80 | opcode, 127, length)

    def _write_locked(self, opcode: int, payload: bytes | bytearray) -> None:
        header = self._frame_header(opcode, len(payload))
        try:
            if len(payload) <= _COALESCE_BYTES:
                # One write, one packet: the socket sends without Nagle's
                # delay, so a header written alone would go out on its own.
                self._writer.write(header + payload)
            else:
                self._writer.write(header)
                self._writer.write(payload)
            self._writer.flush()
        except OSError as error:
            self._mark_lost()
            raise WebSocketClosed(1006, "connection lost") from error

    def _write_frame(self, opcode: int, payload: bytes | bytearray) -> None:
        with self._write_lock:
            if self._closed.is_set():
                raise WebSocketClosed(self.close_code or CLOSE_GOING_AWAY, "connection closed")
            self._write_locked(opcode, payload)

    def send_text(self, text: str) -> None:
        """Send one text message (one frame)."""
        self._write_frame(OPCODE_TEXT, text.encode("utf-8"))

    def send_binary(self, payload: bytes | bytearray) -> None:
        """Send one binary message (one frame), its bytes as they are."""
        self._write_frame(OPCODE_BINARY, payload)

    def ping(self) -> None:
        """Send one ping; the browser answers it with a pong by itself."""
        self._write_frame(OPCODE_PING, b"")

    def close(self, code: int, reason: str = "") -> None:
        """Send a close frame once, then shut the write side down.

        Nothing is sent after the close frame, and the peer reads the end of
        the stream right after it. The read side stays open for the
        connection's own thread to drain (:meth:`finish`): a thread blocked
        reading sees the close within :data:`_CLOSED_CHECK_SECONDS`, even when
        the peer sends nothing more. Safe from any thread, and it never
        blocks longer than the close frame's own write.
        """
        # Cut on a character boundary: a browser fails a connection whose close
        # reason is not valid UTF-8.
        encoded = (
            reason.encode("utf-8")[:_CLOSE_REASON_BYTES].decode("utf-8", "ignore").encode("utf-8")
        )
        with self._write_lock:
            if self._closed.is_set():
                return
            try:
                self._write_locked(OPCODE_CLOSE, struct.pack("!H", code) + encoded)
            except WebSocketClosed:
                return
            self.close_code = code
            self._drain_deadline = time.monotonic() + CLOSE_DRAIN_SECONDS
            self._closed.set()
            with contextlib.suppress(OSError):
                self._connection.shutdown(socket.SHUT_WR)

    def _mark_lost(self) -> None:
        """The connection failed under a read or a write: there is nothing to close."""
        if self._closed.is_set():
            return
        self.close_code = 1006
        self._closed.set()
        with contextlib.suppress(OSError):
            self._connection.shutdown(socket.SHUT_RDWR)

    def finish(self) -> None:
        """Drain the connection after its close; the reading thread's last call.

        What the peer still sends -- its close frame, the rest of a message
        the server refused from its header -- is read and discarded until it
        closes its side, or until :data:`CLOSE_DRAIN_SECONDS` after the close
        frame went out. The socket is then closed with no input unread, so it
        ends with a FIN, not a reset. A connection that was lost has nothing
        to drain.
        """
        self._received.clear()
        deadline = self._drain_deadline
        if deadline is None:
            return
        try:
            while (left := deadline - time.monotonic()) > 0:
                if not self._poller(left) or not self._connection.recv_into(self._scratch):
                    return
        except OSError:
            # Reset by the peer: the close that follows is all that is left.
            return

    # ------------------------------------------------------------------
    # Reading
    # ------------------------------------------------------------------

    def _receive_into(self, target: memoryview) -> int:
        """Receive at least one byte into ``target``; how many.

        Once the connection is closed -- by this thread, or by another while
        this one waited -- nothing more is read as a message: what the peer
        still sends is :meth:`finish`'s to drain.
        """
        while True:
            if self._closed.is_set():
                raise WebSocketClosed(self.close_code or 1006, "connection closed")
            try:
                if not self._poller(_CLOSED_CHECK_SECONDS):
                    continue
                count = self._connection.recv_into(target)
            except (OSError, ValueError) as error:
                self._mark_lost()
                raise WebSocketClosed(1006, "connection lost") from error
            if not count:
                self._mark_lost()
                raise WebSocketClosed(1006, "connection closed mid-frame")
            return count

    def _read_into(self, target: memoryview) -> None:
        filled = min(len(self._received), len(target))
        if filled:
            target[:filled] = self._received[:filled]
            del self._received[:filled]
        while filled < len(target):
            filled += self._receive_into(target[filled:])

    def _read_bytes(self, size: int) -> bytes:
        # Frame headers are a few bytes each: read ahead, so a run of small
        # frames is not a system call per header field.
        while len(self._received) < size:
            count = self._receive_into(self._scratch)
            self._received += self._scratch[:count]
        data = bytes(self._received[:size])
        del self._received[:size]
        return data

    def _read_frame_header(self) -> _Frame:
        first, second = self._read_bytes(2)
        if first & 0x70:
            raise WebSocketProtocolError(CLOSE_PROTOCOL_ERROR, "No extension was negotiated")
        if not second & 0x80:
            raise WebSocketProtocolError(CLOSE_PROTOCOL_ERROR, "Client frames must be masked")
        fin = bool(first & 0x80)
        opcode = first & 0x0F
        length = second & 0x7F
        if length == 126:
            (length,) = struct.unpack("!H", self._read_bytes(2))
        elif length == 127:
            (length,) = struct.unpack("!Q", self._read_bytes(8))
            if length >> 63:
                raise WebSocketProtocolError(CLOSE_PROTOCOL_ERROR, "Frame length is out of range")
        mask = self._read_bytes(4)
        if opcode >= 0x8 and (not fin or length > 125):
            raise WebSocketProtocolError(
                CLOSE_PROTOCOL_ERROR, "Control frames must be whole and at most 125 bytes"
            )
        return _Frame(fin, opcode, length, mask)

    def _control(self, frame: _Frame) -> None:
        payload = bytearray(frame.remaining)
        view = memoryview(payload)
        self._read_into(view)
        unmask_in_place(view, frame.mask, 0)
        if frame.opcode == OPCODE_PING:
            self._write_frame(OPCODE_PONG, bytes(payload))
        elif frame.opcode == OPCODE_CLOSE:
            code = struct.unpack("!H", payload[:2])[0] if len(payload) >= 2 else 1005
            reason = bytes(payload[2:]).decode("utf-8", "replace")
            self.close(code if len(payload) >= 2 else CLOSE_NORMAL)
            raise WebSocketClosed(code, reason)
        elif frame.opcode != OPCODE_PONG:
            raise WebSocketProtocolError(CLOSE_PROTOCOL_ERROR, "Unknown control frame")

    def _next_data_frame(self, *, continuation: bool) -> _Frame:
        while True:
            frame = self._read_frame_header()
            if frame.opcode >= 0x8:
                self._control(frame)
                continue
            if continuation != (frame.opcode == OPCODE_CONTINUATION):
                raise WebSocketProtocolError(
                    CLOSE_PROTOCOL_ERROR,
                    "A fragmented message was interrupted"
                    if continuation
                    else "A continuation frame arrived outside a message",
                )
            if frame.opcode not in (OPCODE_CONTINUATION, OPCODE_TEXT, OPCODE_BINARY):
                raise WebSocketProtocolError(CLOSE_PROTOCOL_ERROR, "Unknown data frame")
            return frame

    def receive(self) -> str | BinaryMessage:
        """The next message: its text, or a :class:`BinaryMessage` to read.

        Raises :class:`WebSocketClosed` when the connection closes (a close
        frame is answered first), or was closed by this side already: nothing
        a peer sends after its close is read as a message.
        Raises :class:`WebSocketProtocolError` when the peer breaks the
        framing rules.
        """
        if self._closed.is_set():
            raise WebSocketClosed(self.close_code or CLOSE_GOING_AWAY, "connection closed")
        frame = self._next_data_frame(continuation=False)
        if frame.opcode == OPCODE_BINARY:
            return BinaryMessage(self, frame)
        parts = bytearray()
        while True:
            if len(parts) + frame.remaining > self._max_text_bytes:
                raise WebSocketProtocolError(
                    CLOSE_MESSAGE_TOO_BIG,
                    f"Text messages are at most {self._max_text_bytes} {plural(self._max_text_bytes, 'byte')}",
                )
            start = len(parts)
            parts.extend(bytes(frame.remaining))
            with memoryview(parts) as whole, whole[start:] as view:
                self._read_into(view)
                unmask_in_place(view, frame.mask, 0)
            if frame.fin:
                break
            frame = self._next_data_frame(continuation=True)
        try:
            return parts.decode("utf-8")
        except UnicodeDecodeError as error:
            raise WebSocketProtocolError(CLOSE_INVALID_DATA, "Text is not UTF-8") from error
