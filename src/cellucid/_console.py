"""Encoding-safe, at-once writes for user-facing console output."""

from __future__ import annotations

import sys
from typing import TextIO


def _text_for_console_stream(text: str, stream: object) -> str:
    """Return text that the stream can encode without changing stream policy."""
    encoding = getattr(stream, "encoding", None)
    if not encoding:
        return text
    try:
        text.encode(encoding, errors="strict")
    except UnicodeEncodeError:
        return text.encode(encoding, errors="backslashreplace").decode(encoding)
    return text


def console_print(
    *values: object,
    sep: str | None = " ",
    end: str | None = "\n",
    file: TextIO | None = None,
) -> None:
    """Print once, escaping only characters unsupported by the target stream,
    and flush it: the line reaches the person at once.

    Python block-buffers standard output when it is not a terminal -- a file,
    a pipe, a process manager, ``nohup``, an IDE's run panel -- so a line
    written and not flushed waits there while a server runs: ``cellucid serve
    data.h5ad > serve.out`` printed its steps but not its Viewer URL until it
    stopped. Every line the person must see goes through here.

    Python's Windows console streams can use single-byte encodings such as
    CP1252. User-controlled names and paths, as well as decorative output
    glyphs, must not turn a successful command into ``UnicodeEncodeError``.
    Backslash escaping is lossless, ASCII-readable, and deterministic without
    mutating the process-wide stream configuration.
    """
    if sep is None:
        sep = " "
    elif not isinstance(sep, str):
        raise TypeError("sep must be None or a string")
    if end is None:
        end = "\n"
    elif not isinstance(end, str):
        raise TypeError("end must be None or a string")

    stream = sys.stdout if file is None else file
    text = sep.join(str(value) for value in values) + end
    stream.write(_text_for_console_stream(text, stream))
    stream.flush()
