"""
Local agent connections: how ``cellucid mcp`` reaches a Cellucid viewer.

An agent app on this computer (Claude Code, Claude Desktop, Cursor) runs
``cellucid mcp`` (:mod:`cellucid.mcp`), a local MCP server over stdio. Each of
its tool calls is carried to a viewer page over that viewer's channel as an
agent command (``agent_list``, ``agent_call``, ``agent_cancel``,
``agent_answer``; ``viewer_protocol.json``), and the page runs it exactly as
it runs a connector agent's call: as the agent, through a connection of kind
``local``, under the person's permissions, asked in the confirm card, one
turn the person can Stop and undo. The page answers with agent events
(``agent_tools``, ``agent_result``, ``agent_asking``, ``agent_answered``).

Two ways to a viewer, one interface (:class:`ViewerLink`):

- **In the same process** (``cellucid mcp --open <file>`` serves the viewer
  itself): the viewer's :class:`ViewerAgentPort` sends the commands on its
  channel and hands each answer to whoever waits for it.
- **In another process** (a notebook, a script, an R session): every viewer
  with a channel is published to agents on this computer. Its process runs one
  loopback HTTP endpoint (:class:`_AgentEndpoint`, started with the first
  viewer, on a port the system picks) and writes one record per viewer into
  the agents folder (:func:`agents_dir`): where the endpoint is and the secret
  token it takes, readable by this user only. ``cellucid mcp`` reads the
  records (:func:`discover_viewers`) and drives one through
  :class:`RemoteViewerLink`. cellucid-r publishes its viewers the same way,
  so either package's ``cellucid mcp`` drives either package's viewer.

The endpoint (``<endpoint>`` is ``http://127.0.0.1:<port>/_cellucid/agent/<viewer id>``):

- ``GET <endpoint>``: the viewer (:meth:`ViewerAgentPort.describe`).
- ``POST <endpoint>/command``: one agent command, exactly as the viewer
  protocol has it; ``202`` once sent, ``409`` with the reason when no page is
  connected to the viewer, ``400`` for a command that is not one.
- ``GET <endpoint>/events?consumer=<connection id>&wait=<seconds>``: the agent
  events for that consumer (the answers to its commands), oldest first, as
  ``{"events": [...]}``; it waits up
  to ``wait`` seconds (at most :data:`EVENT_WAIT_LIMIT_SECONDS`) for one. A
  command the page refused is ``{"type": "agent_failed", "requestId",
  "reason"}``, the viewer stopping ``{"type": "viewer_stopped", "reason"}``.

Every request carries ``Authorization: Bearer <token>``; one from a browser
(any ``Origin`` header), from another host, or naming another ``Host`` is
refused. The token is compared in constant time.
"""

from __future__ import annotations

import contextlib
import hmac
import json
import logging
import os
import secrets
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import deque
from collections.abc import Callable
from datetime import UTC, datetime
from http.server import BaseHTTPRequestHandler
from pathlib import Path
from typing import Any, Protocol

from ._server_base import CellucidThreadingHTTPServer

logger = logging.getLogger("cellucid.agents")

#: The version of a viewer record and of the endpoint's answers.
AGENT_RECORD_SCHEMA = 2

#: The route every viewer's endpoint is under.
AGENT_ROUTE = "_cellucid/agent"

#: The longest a ``GET .../events`` waits for an event before answering none.
EVENT_WAIT_LIMIT_SECONDS = 25.0

#: How long a consumer that stopped asking for its events keeps them: its
#: ``cellucid mcp`` ended without saying so (killed, the app closed).
CONSUMER_IDLE_SECONDS = 120.0

#: The largest command body the endpoint reads: one tool call written by a
#: model (a model writes at most a few hundred kilobytes in one answer).
COMMAND_BODY_LIMIT_BYTES = 4 * 1024 * 1024

#: How long a connection to the endpoint may sit idle before its thread
#: closes it: longer than a client takes between one request and its next
#: (``cellucid mcp`` asks again at once), as the data server's idle bound.
CONNECTION_IDLE_SECONDS = 120.0

#: What a call's listener hears, in order: ``asking`` and ``answered`` while
#: the person is asked (any number of times), then exactly one of ``result``
#: (the MCP tool result), ``failed`` (the page refused the command:
#: ``reason``) or ``stopped`` (the viewer stopped or its page went: ``reason``).
CALL_EVENT_KINDS = ("asking", "answered", "result", "failed", "stopped")

CallListener = Callable[[str, dict[str, Any]], None]


class ViewerUnavailableError(RuntimeError):
    """No page shows the viewer now, or the viewer stopped: one sentence that says which."""


class ViewerLink(Protocol):
    """One viewer an agent connection drives, wherever it runs."""

    def describe(self) -> dict[str, Any]: ...

    def connected(self) -> bool: ...

    def list_tools(self, connection: dict[str, Any], timeout: float) -> list[dict[str, Any]]: ...

    def start_call(
        self,
        request_id: str,
        tool: str,
        args: dict[str, Any],
        agent: str,
        connection: dict[str, Any],
        listener: CallListener,
    ) -> None: ...

    def cancel(self, request_id: str) -> None: ...

    def answer(self, request_id: str, allowed: bool) -> None: ...

    def close(self) -> None: ...


#: A call whose viewer page went (reloaded or closed) before it was answered.
PAGE_GONE = (
    "The Cellucid viewer's page went (it was reloaded or closed) before this call ended, so whether it "
    "was done is not known. Read the history (history_get) before trying again."
)

#: A call the viewer page refused to read (an older page): ``{reason}`` is the page's.
NOT_TAKEN = "The Cellucid viewer page did not take the call: {reason} Ask the person to reload the viewer's page."

#: A call whose viewer stopped before it was answered.
STOPPED = "The Cellucid viewer stopped ({reason}); nothing more was done."

#: The agents folder shared with others: refused, never used.
FOLDER_UNSAFE = (
    "The agents folder {folder} is not this user's alone ({why}), so no viewer is published there or "
    "read from it. Make it yours and not writable by others (chmod 700), or set CELLUCID_AGENTS_DIR."
)

#: Set to ``off`` (or ``0``, ``false``, ``no``), viewers are not published to agents on this computer.
OPT_OUT_VARIABLE = "CELLUCID_AGENTS"


def agents_disabled() -> bool:
    """Whether this user turned publication off (``CELLUCID_AGENTS=off``)."""
    return os.environ.get(OPT_OUT_VARIABLE, "").strip().lower() in {"off", "0", "false", "no"}


def not_open_sentence(description: dict[str, Any]) -> str:
    """What an agent reads when no page shows a viewer now: open its page
    again (the address its program printed: it carries the viewer's secret,
    which an agent is never given), or show its notebook cell again."""
    name = (description.get("dataset") or {}).get("name") or "the dataset"
    url = description.get("url")
    if url:
        return (
            f"No page shows the Cellucid viewer of {name} (served at {url}) now. Ask the person to open "
            "its page again (the viewer address its program printed, or reload the page), then try again."
        )
    return (
        f"No page shows the Cellucid viewer of {name} now. Ask the person to show it again "
        "(run the notebook cell that shows it), then try again."
    )


# =============================================================================
# One viewer's end, in its own process
# =============================================================================


class ViewerAgentPort:
    """A viewer's end of local agent connections (:class:`ViewerLink`).

    ``send`` sends one command on the viewer's channel (it is checked and
    serialized there); ``connected`` says whether a page holds the channel
    now; ``describe`` names the viewer (:meth:`describe`). The viewer hands
    every agent event to :meth:`deliver`, a ``command_error`` of an agent
    command to :meth:`command_failed`, and calls :meth:`stopped` when it
    stops; each reaches whoever waits for it, on the channel's thread.
    """

    def __init__(
        self,
        *,
        send: Callable[[dict[str, Any]], object],
        connected: Callable[[], bool],
        describe: Callable[[], dict[str, Any]],
    ) -> None:
        self._send = send
        self._connected = connected
        self._describe = describe
        self._lock = threading.Lock()
        self._calls: dict[str, CallListener] = {}
        self._lists: dict[str, list[Any]] = {}
        self._stopped: str | None = None

    # -- the viewer's side -------------------------------------------------

    def deliver(self, event: str, data: dict[str, Any]) -> None:
        """One agent event from the page (``data`` without its ``type``)."""
        request_id = data["requestId"]
        if event == "agent_tools":
            self._settle_list(request_id, ("tools", data["tools"]))
            return
        with self._lock:
            listener = self._calls.get(request_id)
            if event == "agent_result":
                self._calls.pop(request_id, None)
        if listener is None:
            return
        if event == "agent_result":
            self._hear(listener, "result", {"result": data["result"]})
        elif event == "agent_asking":
            self._hear(listener, "asking", {"message": data["message"], "detail": data["detail"]})
        elif event == "agent_answered":
            self._hear(listener, "answered", {"allowed": data["allowed"]})

    def command_failed(self, data: dict[str, Any]) -> None:
        """The page refused an agent command (its ``command_error``): an old
        page, or a command it could not read."""
        request_id = data.get("requestId")
        if request_id is None:
            return
        reason = NOT_TAKEN.format(reason=data["reason"])
        if data["command"] == "agent_list":
            self._settle_list(request_id, ("failed", reason))
            return
        with self._lock:
            listener = self._calls.pop(request_id, None)
        if listener is not None:
            self._hear(listener, "failed", {"reason": reason})

    def page_gone(self) -> None:
        """The page holding the channel went (closed, reloaded, its link
        dropped): no answer of a call it was running will come, so each
        fails now with the sentence an agent reads (:data:`PAGE_GONE`)."""
        with self._lock:
            calls = list(self._calls.values())
            self._calls.clear()
            lists = list(self._lists)
        for listener in calls:
            self._hear(listener, "failed", {"reason": PAGE_GONE})
        for request_id in lists:
            self._settle_list(request_id, ("failed", PAGE_GONE))

    def stopped(self, reason: str) -> None:
        """The viewer stopped: everything waiting ends with ``reason``."""
        with self._lock:
            self._stopped = reason
            calls = list(self._calls.values())
            self._calls.clear()
            lists = list(self._lists)
        for listener in calls:
            self._hear(listener, "stopped", {"reason": reason})
        for request_id in lists:
            self._settle_list(request_id, ("stopped", reason))

    # -- the agent's side (ViewerLink) ---------------------------------------

    def describe(self) -> dict[str, Any]:
        return {**self._describe(), "connected": self.connected()}

    def connected(self) -> bool:
        if self._stopped is not None:
            return False
        try:
            return bool(self._connected())
        except LookupError:
            return False

    def _require_open(self) -> None:
        if self._stopped is not None:
            raise ViewerUnavailableError(f"The Cellucid viewer stopped ({self._stopped}).")
        if not self.connected():
            raise ViewerUnavailableError(not_open_sentence(self._describe()))

    def list_tools(self, connection: dict[str, Any], timeout: float) -> list[dict[str, Any]]:
        """The tools as ``connection`` may use them, from the page."""
        self._require_open()
        request_id = f"{connection['id']}-list-{secrets.token_hex(4)}"
        slot: list[Any] = [threading.Event(), None]
        with self._lock:
            self._lists[request_id] = slot
        try:
            self._send({"type": "agent_list", "requestId": request_id, "connection": dict(connection)})
            if not slot[0].wait(timeout):
                raise TimeoutError("The Cellucid viewer did not list its tools in time.")
        finally:
            with self._lock:
                self._lists.pop(request_id, None)
        kind, value = slot[1]
        if kind == "tools":
            return list(value)
        if kind == "failed":
            raise ViewerUnavailableError(f"The Cellucid viewer page did not list its tools: {value}")
        raise ViewerUnavailableError(f"The Cellucid viewer stopped ({value}).")

    def start_call(
        self,
        request_id: str,
        tool: str,
        args: dict[str, Any],
        agent: str,
        connection: dict[str, Any],
        listener: CallListener,
    ) -> None:
        """Send one call; ``listener`` hears its events (:data:`CALL_EVENT_KINDS`)."""
        self._require_open()
        with self._lock:
            self._calls[request_id] = listener
        try:
            self._send(
                {
                    "type": "agent_call",
                    "requestId": request_id,
                    "tool": tool,
                    "args": args,
                    "agent": agent,
                    "connection": dict(connection),
                }
            )
        except BaseException:
            with self._lock:
                self._calls.pop(request_id, None)
            raise

    def cancel(self, request_id: str) -> None:
        with self._lock:
            waiting = self._calls.pop(request_id, None) is not None
        if waiting and self.connected():
            self._send({"type": "agent_cancel", "requestId": request_id})

    def answer(self, request_id: str, allowed: bool) -> None:
        if self.connected():
            self._send({"type": "agent_answer", "requestId": request_id, "allowed": bool(allowed)})

    def close(self) -> None:
        """Nothing to close: the viewer owns its channel."""

    # -- helpers -------------------------------------------------------------

    def _settle_list(self, request_id: str, outcome: tuple[str, Any]) -> None:
        with self._lock:
            slot = self._lists.get(request_id)
        if slot is not None and slot[1] is None:
            slot[1] = outcome
            slot[0].set()

    @staticmethod
    def _hear(listener: CallListener, kind: str, data: dict[str, Any]) -> None:
        try:
            listener(kind, data)
        except Exception:
            logger.exception("An agent call listener failed")


# =============================================================================
# The agents folder: one record per published viewer
# =============================================================================


def agents_dir() -> Path:
    """The folder of this user's published viewers.

    ``$CELLUCID_AGENTS_DIR`` when set (an absolute path; a relative one is
    refused with :class:`AgentsFolderError`), else
    ``<user cache>/cellucid/agents``: ``%LOCALAPPDATA%`` on Windows,
    ``~/Library/Caches`` on macOS, ``$XDG_CACHE_HOME`` or ``~/.cache``
    elsewhere -- the same folder cellucid-r reads and writes, so either
    package's ``cellucid mcp`` finds either package's viewers.
    """
    override = os.environ.get("CELLUCID_AGENTS_DIR", "")
    if override:
        path = Path(override)
        if not path.is_absolute():
            # Refused as a folder shared with others is: nothing is published
            # or read, and this is the sentence an agent reads (as cellucid-r's).
            raise AgentsFolderError("CELLUCID_AGENTS_DIR must be an absolute path")
        return path
    if sys.platform == "win32":
        local = os.environ.get("LOCALAPPDATA", "")
        base = Path(local) if local and Path(local).is_absolute() else Path.home() / "AppData" / "Local"
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        xdg = os.environ.get("XDG_CACHE_HOME", "")
        base = Path(xdg) if xdg and Path(xdg).is_absolute() else Path.home() / ".cache"
    return base / "cellucid" / "agents"


#: The keys of a viewer record, exactly (cellucid-r writes the same).
RECORD_FIELDS = frozenset(
    {"schema", "viewer", "runtime", "host", "pid", "started", "dataset", "url", "endpoint", "token"}
)

#: The hosts an endpoint may name: this computer's loopback, nothing else.
LOOPBACK_HOSTS = frozenset({"127.0.0.1", "localhost", "::1"})


class AgentsFolderError(RuntimeError):
    """The agents folder is refused (not this user's alone, or named by a
    relative ``CELLUCID_AGENTS_DIR``): one sentence that says why."""


def _this_host() -> str:
    import socket

    return socket.gethostname()


def _require_private_folder(folder: Path) -> None:
    """On POSIX, the agents folder must be this user's and writable by no one
    else: a record names an endpoint and its secret, and a folder others can
    write could hand ``cellucid mcp`` a record of theirs. A home folder
    shared between machines is fine (records name their host)."""
    if sys.platform == "win32":
        return
    status = folder.stat()
    if status.st_uid != os.getuid():
        raise AgentsFolderError(FOLDER_UNSAFE.format(folder=folder, why="another user owns it"))
    if status.st_mode & 0o022:
        raise AgentsFolderError(FOLDER_UNSAFE.format(folder=folder, why="others may write in it"))


def _private_folder() -> Path:
    """The agents folder, made (0700) when missing, and checked."""
    folder = agents_dir()
    if not folder.exists():
        folder.mkdir(parents=True, exist_ok=True)
        if sys.platform != "win32":
            os.chmod(folder, 0o700)
    _require_private_folder(folder)
    return folder


def _write_record(record: dict[str, Any]) -> Path:
    folder = _private_folder()
    target = folder / f"{record['viewer']}.json"
    descriptor, temporary = tempfile.mkstemp(dir=folder, prefix=".", suffix=".tmp")
    try:
        if sys.platform != "win32":
            os.fchmod(descriptor, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(record, handle, ensure_ascii=False, sort_keys=True)
        os.replace(temporary, target)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(temporary)
        raise
    return target


def _is_loopback_endpoint(endpoint: str) -> bool:
    try:
        parts = urllib.parse.urlsplit(endpoint)
    except ValueError:
        return False
    return parts.scheme == "http" and parts.hostname in LOOPBACK_HOSTS and parts.port is not None


def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    if sys.platform == "win32":
        import ctypes

        kernel32 = ctypes.windll.kernel32  # type: ignore[attr-defined]
        handle = kernel32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if not handle:
            return False
        kernel32.CloseHandle(handle)
        return True
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _text(value: object) -> bool:
    return type(value) is str and value != ""


def _read_record(path: Path) -> dict[str, Any] | None:
    """A viewer record, or None for anything that is not exactly one: every
    field of its type, its endpoint on this computer's loopback."""
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if type(record) is not dict or set(record) != RECORD_FIELDS or record["schema"] != AGENT_RECORD_SCHEMA:
        return None
    dataset = record["dataset"]
    if not (
        all(_text(record[key]) for key in ("viewer", "runtime", "host", "started", "endpoint", "token"))
        and type(record["pid"]) is int
        and (record["url"] is None or _text(record["url"]))
        and (dataset is None or (type(dataset) is dict and set(dataset) == {"id", "name"} and all(_text(dataset[key]) for key in ("id", "name"))))
        and _is_loopback_endpoint(record["endpoint"])
    ):
        return None
    return record


def discover_viewers() -> list[dict[str, Any]]:
    """Every published viewer of this user on this computer whose process
    still runs, newest first. Records of this computer's processes that ended
    are removed; another computer's (a home folder shared between machines)
    are left alone. A folder that is not this user's alone is refused
    (:class:`AgentsFolderError`)."""
    folder = agents_dir()
    if not folder.exists():
        return []
    _require_private_folder(folder)
    try:
        paths = sorted(folder.glob("*.json"))
    except OSError:
        return []
    host = _this_host()
    records = []
    for path in paths:
        record = _read_record(path)
        if record is None or record["host"] != host:
            continue
        if not _pid_alive(record["pid"]):
            with contextlib.suppress(OSError):
                path.unlink()
            continue
        records.append(record)
    records.sort(key=lambda record: record["started"], reverse=True)
    return records


# =============================================================================
# The process's endpoint: its viewers, reachable by agents on this computer
# =============================================================================


class _Consumer:
    __slots__ = ("events", "requests", "seen")

    def __init__(self) -> None:
        self.events: deque[dict[str, Any]] = deque()
        self.requests: set[str] = set()
        self.seen = time.monotonic()


class _PublishedViewer:
    """One viewer behind the endpoint: its port, its record, and the queue of
    each consumer (one ``cellucid mcp`` connection) of its events."""

    def __init__(self, viewer_id: str, port: ViewerAgentPort, runtime: str) -> None:
        self.viewer_id = viewer_id
        self.port = port
        self.runtime = runtime
        self.token = secrets.token_hex(32)
        self.lock = threading.Condition()
        self.consumers: dict[str, _Consumer] = {}
        self.record_path: Path | None = None
        self.stopped: str | None = None

    def sweep(self, keep: str | None = None) -> None:
        """Drop the consumers silent for :data:`CONSUMER_IDLE_SECONDS` (their
        ``cellucid mcp`` ended without saying so): their calls are stopped in
        the page (a question withdrawn, a change rolled back) and what waited
        for them is let go. The caller holds ``lock``."""
        now = time.monotonic()
        for other, stale in list(self.consumers.items()):
            if now - stale.seen > CONSUMER_IDLE_SECONDS and other != keep:
                for request_id in stale.requests:
                    with contextlib.suppress(Exception):
                        self.port.cancel(request_id)
                del self.consumers[other]

    def consumer(self, name: str) -> _Consumer:
        now = time.monotonic()
        self.sweep(keep=name)
        entry = self.consumers.get(name)
        if entry is None:
            entry = self.consumers[name] = _Consumer()
        entry.seen = now
        return entry

    def push(self, name: str, event: dict[str, Any]) -> None:
        with self.lock:
            entry = self.consumers.get(name)
            if entry is None:
                return
            entry.events.append(event)
            if event["type"] in {"agent_result", "agent_failed"}:
                entry.requests.discard(event["requestId"])
            self.lock.notify_all()

    def stop(self, reason: str) -> None:
        with self.lock:
            self.stopped = reason
            for entry in self.consumers.values():
                entry.events.append({"type": "viewer_stopped", "reason": reason})
            self.lock.notify_all()


class _AgentEndpoint:
    """The process's one loopback endpoint for agents (started with its first published viewer)."""

    def __init__(self) -> None:
        self.viewers: dict[str, _PublishedViewer] = {}
        self.lock = threading.Lock()
        endpoint = self

        class Handler(_EndpointHandler):
            owner = endpoint

        # The package's server: a client that goes mid-request is no error
        # printed in a notebook, and on Windows its port is bound exclusively.
        self.server = CellucidThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.port = int(self.server.server_address[1])
        self.thread = threading.Thread(target=self.server.serve_forever, name="cellucid-agents", daemon=True)
        self.thread.start()
        # Consumers that went silent are swept on a clock, not only when another one asks.
        self.sweeper = threading.Thread(target=self._sweep_forever, name="cellucid-agents-sweep", daemon=True)
        self.sweeper.start()

    def _sweep_forever(self) -> None:
        while True:
            time.sleep(CONSUMER_IDLE_SECONDS / 4)
            with self.lock:
                viewers = list(self.viewers.values())
            for published in viewers:
                with published.lock:
                    published.sweep()

    def url_of(self, viewer_id: str) -> str:
        return f"http://127.0.0.1:{self.port}/{AGENT_ROUTE}/{viewer_id}"


_endpoint: _AgentEndpoint | None = None
_endpoint_lock = threading.Lock()


def publish_viewer(viewer_id: str, port: ViewerAgentPort, *, runtime: str) -> None:
    """Publish a viewer to agents on this computer: its record in the agents
    folder, its commands and events on this process's endpoint. A folder that
    cannot be written leaves the viewer unpublished, said in the log: the
    viewer works on, and ``cellucid mcp --open`` still reaches its own.
    ``CELLUCID_AGENTS=off`` publishes nothing."""
    global _endpoint
    if agents_disabled():
        return
    try:
        with _endpoint_lock:
            if _endpoint is None:
                _endpoint = _AgentEndpoint()
            endpoint = _endpoint
        published = _PublishedViewer(viewer_id, port, runtime)
        with endpoint.lock:
            endpoint.viewers[viewer_id] = published
        description = port.describe()
        published.record_path = _write_record(
            {
                "schema": AGENT_RECORD_SCHEMA,
                "viewer": viewer_id,
                "runtime": runtime,
                "host": _this_host(),
                "pid": os.getpid(),
                "started": datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%S.%fZ"),
                "dataset": description.get("dataset"),
                "url": description.get("url"),
                "endpoint": endpoint.url_of(viewer_id),
                "token": published.token,
            }
        )
    except (OSError, ValueError, AgentsFolderError) as error:
        # Not published: nothing of it stays behind the endpoint.
        endpoint_now = _endpoint
        if endpoint_now is not None:
            with endpoint_now.lock:
                endpoint_now.viewers.pop(viewer_id, None)
        logger.warning("This viewer is not published to agents on this computer: %s", error)
    except Exception:
        # Publication is the agents' way in, never the viewer's own work: a
        # defect here is reported in full and the viewer goes on.
        logger.warning("This viewer could not be published to agents on this computer", exc_info=True)


def withdraw_viewer(viewer_id: str, reason: str = "the viewer was stopped") -> None:
    """Withdraw a published viewer: its record goes, its consumers hear it stopped."""
    endpoint = _endpoint
    if endpoint is None:
        return
    with endpoint.lock:
        published = endpoint.viewers.pop(viewer_id, None)
    if published is None:
        return
    published.stop(reason)
    if published.record_path is not None:
        with contextlib.suppress(OSError):
            published.record_path.unlink()


class _EndpointHandler(BaseHTTPRequestHandler):
    owner: _AgentEndpoint
    server_version = "cellucid-agents"
    protocol_version = "HTTP/1.1"
    # A connection that sends nothing (its next request) for this long is
    # closed, so it holds no thread: anyone on this computer can open one
    # without the token. A waiting GET .../events does not read the socket,
    # so the wait is not cut short by this.
    timeout = CONNECTION_IDLE_SECONDS

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002 - the base class's name
        logger.debug("agents endpoint: " + format, *args)

    def _answer(self, status: int, value: dict[str, Any]) -> None:
        body = json.dumps(value, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _viewer(self) -> tuple[_PublishedViewer, str, dict[str, list[str]]] | None:
        """The viewer a request names, after its checks; None once answered."""
        if self.headers.get("Origin") is not None:
            self._answer(403, {"error": "Browsers cannot use this endpoint."})
            return None
        host = (self.headers.get("Host") or "").rsplit(":", 1)[0]
        if host not in {"127.0.0.1", "localhost"} or self.client_address[0] != "127.0.0.1":
            self._answer(403, {"error": "Only this computer can use this endpoint."})
            return None
        parts = urllib.parse.urlsplit(self.path)
        segments = parts.path.strip("/").split("/")
        if len(segments) < 3 or "/".join(segments[:2]) != AGENT_ROUTE:
            self._answer(404, {"error": "Not found."})
            return None
        with self.owner.lock:
            published = self.owner.viewers.get(segments[2])
        if published is None:
            self._answer(404, {"error": "This viewer is not published (it stopped)."})
            return None
        given = (self.headers.get("Authorization") or "").removeprefix("Bearer ")
        if not hmac.compare_digest(given.encode("utf-8"), published.token.encode("utf-8")):
            self._answer(401, {"error": "The token does not match this viewer."})
            return None
        return published, "/".join(segments[3:]), urllib.parse.parse_qs(parts.query)

    def do_GET(self) -> None:  # noqa: N802 - the base class's name
        found = self._viewer()
        if found is None:
            return
        published, rest, query = found
        if rest == "":
            self._answer(200, {"schema": AGENT_RECORD_SCHEMA, "viewer": published.port.describe()})
            return
        if rest != "events":
            self._answer(404, {"error": "Not found."})
            return
        consumer = (query.get("consumer") or [""])[0]
        try:
            wait = min(max(float((query.get("wait") or ["0"])[0]), 0.0), EVENT_WAIT_LIMIT_SECONDS)
        except ValueError:
            wait = 0.0
        if not consumer:
            self._answer(400, {"error": "Name the consumer (?consumer=<connection id>)."})
            return
        deadline = time.monotonic() + wait
        with published.lock:
            entry = published.consumer(consumer)
            while not entry.events and published.stopped is None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    break
                published.lock.wait(remaining)
            events = list(entry.events)
            entry.events.clear()
            if published.stopped is not None and not any(event["type"] == "viewer_stopped" for event in events):
                events.append({"type": "viewer_stopped", "reason": published.stopped})
        self._answer(200, {"events": events})

    def do_POST(self) -> None:  # noqa: N802 - the base class's name
        found = self._viewer()
        if found is None:
            return
        published, rest, _query = found
        if rest != "command":
            self._answer(404, {"error": "Not found."})
            return
        try:
            length = int(self.headers.get("Content-Length") or "0")
        except ValueError:
            length = 0
        if length <= 0 or length > COMMAND_BODY_LIMIT_BYTES:
            self._answer(400, {"error": "A command is one JSON object."})
            return
        try:
            command = json.loads(self.rfile.read(length).decode("utf-8"))
        except (UnicodeDecodeError, ValueError, RecursionError):
            self._answer(400, {"error": "A command is one JSON object."})
            return
        try:
            self._run(published, command)
        except ViewerUnavailableError as error:
            self._answer(409, {"error": str(error)})
            return
        except (TypeError, ValueError, KeyError) as error:
            self._answer(400, {"error": str(error)})
            return
        except RecursionError:
            self._answer(400, {"error": "The command is nested more deeply than this package reads."})
            return
        self._answer(202, {})

    @staticmethod
    def _run(published: _PublishedViewer, command: Any) -> None:
        if type(command) is not dict or type(command.get("type")) is not str:
            raise TypeError("A command is one JSON object with its type.")
        kind = command["type"]
        request_id = command.get("requestId")
        if type(request_id) is not str:
            raise TypeError("A command names its requestId.")
        port = published.port
        if kind == "agent_list":
            connection = command["connection"]
            name = connection["id"]
            with published.lock:
                published.consumer(name)

            def list_in_background() -> None:
                try:
                    tools = port.list_tools(connection, timeout=60.0)
                except (ViewerUnavailableError, TimeoutError) as error:
                    published.push(name, {"type": "agent_failed", "requestId": request_id, "reason": str(error)})
                    return
                published.push(name, {"type": "agent_tools", "requestId": request_id, "tools": tools})

            threading.Thread(target=list_in_background, name="cellucid-agent-list", daemon=True).start()
            return
        if kind == "agent_call":
            connection = command["connection"]
            name = connection["id"]
            with published.lock:
                published.consumer(name).requests.add(request_id)

            def listener(event: str, data: dict[str, Any]) -> None:
                if event == "result":
                    published.push(name, {"type": "agent_result", "requestId": request_id, "result": data["result"]})
                elif event == "asking":
                    published.push(name, {"type": "agent_asking", "requestId": request_id, **data})
                elif event == "answered":
                    published.push(name, {"type": "agent_answered", "requestId": request_id, **data})
                elif event == "stopped":
                    published.push(name, {"type": "viewer_stopped", "reason": data["reason"]})
                else:
                    published.push(name, {"type": "agent_failed", "requestId": request_id, "reason": data["reason"]})

            try:
                port.start_call(request_id, command["tool"], command["args"], command["agent"], connection, listener)
            except BaseException:
                with published.lock:
                    published.consumer(name).requests.discard(request_id)
                raise
            return
        if kind == "agent_cancel":
            port.cancel(request_id)
            return
        if kind == "agent_answer":
            port.answer(request_id, command["allowed"] is True)
            return
        raise ValueError(f"{kind!r} is no agent command.")


# =============================================================================
# A viewer in another process
# =============================================================================


_NO_PROXY_OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))


def _endpoint_request(
    record: dict[str, Any], method: str, path: str, body: dict[str, Any] | None = None, timeout: float = 10.0
) -> tuple[int, dict[str, Any]]:
    """One request to a published viewer's endpoint: its status and JSON answer
    (never through a proxy: the endpoint is on this computer)."""
    data = None if body is None else json.dumps(body).encode("utf-8")
    request = urllib.request.Request(
        f"{record['endpoint']}{path}",
        data=data,
        method=method,
        headers={
            "Authorization": f"Bearer {record['token']}",
            **({"Content-Type": "application/json"} if data is not None else {}),
        },
    )
    try:
        with _NO_PROXY_OPENER.open(request, timeout=timeout) as response:
            payload = json.loads(response.read().decode("utf-8") or "{}")
            return response.status, payload if type(payload) is dict else {}
    except urllib.error.HTTPError as error:
        try:
            payload = json.loads(error.read().decode("utf-8") or "{}")
        except ValueError:
            payload = {}
        return error.code, payload if type(payload) is dict else {}


def describe_record(record: dict[str, Any]) -> dict[str, Any]:
    """A published viewer as its process describes it now; from its record,
    and not connected, when its process does not answer."""
    try:
        status, payload = _endpoint_request(record, "GET", "", timeout=5.0)
    except (OSError, ValueError):
        status, payload = 0, {}
    if status == 200 and type(payload.get("viewer")) is dict:
        return dict(payload["viewer"])
    return {
        "viewer": record["viewer"],
        "runtime": record["runtime"],
        "dataset": record["dataset"],
        "url": record["url"],
        "connected": False,
    }


class RemoteViewerLink:
    """A viewer another process published (:class:`ViewerLink`), driven over
    its endpoint. One thread reads the events of this connection (the
    consumer) and hands each to whoever waits for it."""

    def __init__(self, record: dict[str, Any], consumer: str) -> None:
        self.record = record
        self._consumer = consumer
        self._lock = threading.Lock()
        self._calls: dict[str, CallListener] = {}
        self._lists: dict[str, list[Any]] = {}
        self._stopped: str | None = None
        self._closed = threading.Event()
        self._reader = threading.Thread(target=self._read_events, name="cellucid-agent-events", daemon=True)
        self._reader.start()

    def _request(
        self, method: str, path: str, body: dict[str, Any] | None = None, timeout: float = 10.0
    ) -> tuple[int, dict[str, Any]]:
        return _endpoint_request(self.record, method, path, body, timeout)

    def describe(self) -> dict[str, Any]:
        description = describe_record(self.record)
        if self._stopped is not None:
            description["connected"] = False
        return description

    def connected(self) -> bool:
        return self._stopped is None and bool(self.describe().get("connected"))

    def _command(self, command: dict[str, Any]) -> None:
        if self._stopped is not None:
            raise ViewerUnavailableError(f"The Cellucid viewer stopped ({self._stopped}).")
        try:
            status, payload = self._request("POST", "/command", command)
        except (OSError, ValueError) as error:
            raise ViewerUnavailableError(f"The Cellucid viewer's process cannot be reached ({error}).") from error
        if status == 202:
            return
        reason = payload.get("error") or f"its process answered {status}"
        if status == 409:
            raise ViewerUnavailableError(str(reason))
        if status in {401, 404}:
            raise ViewerUnavailableError(f"The Cellucid viewer is gone ({reason}).")
        raise ValueError(str(reason))

    def list_tools(self, connection: dict[str, Any], timeout: float) -> list[dict[str, Any]]:
        request_id = f"{connection['id']}-list-{secrets.token_hex(4)}"
        slot: list[Any] = [threading.Event(), None]
        with self._lock:
            self._lists[request_id] = slot
        try:
            self._command({"type": "agent_list", "requestId": request_id, "connection": dict(connection)})
            if not slot[0].wait(timeout):
                raise TimeoutError("The Cellucid viewer did not list its tools in time.")
        finally:
            with self._lock:
                self._lists.pop(request_id, None)
        kind, value = slot[1]
        if kind == "tools":
            return list(value)
        raise ViewerUnavailableError(str(value))

    def start_call(
        self,
        request_id: str,
        tool: str,
        args: dict[str, Any],
        agent: str,
        connection: dict[str, Any],
        listener: CallListener,
    ) -> None:
        with self._lock:
            self._calls[request_id] = listener
        try:
            self._command(
                {
                    "type": "agent_call",
                    "requestId": request_id,
                    "tool": tool,
                    "args": args,
                    "agent": agent,
                    "connection": dict(connection),
                }
            )
        except BaseException:
            with self._lock:
                self._calls.pop(request_id, None)
            raise

    def cancel(self, request_id: str) -> None:
        with self._lock:
            waiting = self._calls.pop(request_id, None) is not None
        if waiting:
            with contextlib.suppress(ViewerUnavailableError, ValueError):
                self._command({"type": "agent_cancel", "requestId": request_id})

    def answer(self, request_id: str, allowed: bool) -> None:
        with contextlib.suppress(ViewerUnavailableError, ValueError):
            self._command({"type": "agent_answer", "requestId": request_id, "allowed": bool(allowed)})

    def close(self) -> None:
        """Stop reading this viewer's events. A call still waiting can get no
        answer through this link now: it fails with the page-gone sentence
        (whether it was done is not known)."""
        self._closed.set()
        with self._lock:
            calls = list(self._calls.values())
            self._calls.clear()
        for listener in calls:
            ViewerAgentPort._hear(listener, "failed", {"reason": PAGE_GONE})

    def _read_events(self) -> None:
        query = urllib.parse.urlencode({"consumer": self._consumer, "wait": EVENT_WAIT_LIMIT_SECONDS})
        failures = 0
        while not self._closed.is_set():
            try:
                status, payload = self._request("GET", f"/events?{query}", timeout=EVENT_WAIT_LIMIT_SECONDS + 10)
            except (OSError, ValueError):
                status, payload = 0, {}
            if status == 200 and type(payload.get("events")) is list:
                failures = 0
                for event in payload["events"]:
                    if type(event) is dict:
                        self._route(event)
                continue
            if status in {401, 404}:
                self._end("the viewer was stopped")
                return
            failures += 1
            if not _pid_alive(self.record["pid"]):
                self._end("its process ended")
                return
            # Its record went: the viewer was stopped (an R viewer's endpoint
            # goes with its server, so nothing answers to say so).
            if not (agents_dir() / f"{self.record['viewer']}.json").exists():
                self._end("the viewer was stopped")
                return
            # Back off while its process is busy or restarting its endpoint.
            self._closed.wait(min(5.0, 0.2 * failures))

    def _route(self, event: dict[str, Any]) -> None:
        kind = event.get("type")
        if kind == "viewer_stopped":
            self._end(str(event.get("reason") or "the viewer was stopped"))
            return
        request_id = event.get("requestId")
        if type(request_id) is not str:
            return
        with self._lock:
            slot = self._lists.get(request_id)
            listener = self._calls.get(request_id)
            if kind in {"agent_result", "agent_failed"}:
                self._calls.pop(request_id, None)
        if slot is not None:
            if kind == "agent_tools":
                slot[1] = ("tools", event.get("tools") or [])
            else:
                slot[1] = ("failed", event.get("reason") or "the viewer did not list its tools")
            slot[0].set()
            return
        if listener is None:
            return
        if kind == "agent_result":
            ViewerAgentPort._hear(listener, "result", {"result": event.get("result")})
        elif kind == "agent_failed":
            ViewerAgentPort._hear(listener, "failed", {"reason": str(event.get("reason"))})
        elif kind == "agent_asking":
            ViewerAgentPort._hear(listener, "asking", {"message": event.get("message"), "detail": event.get("detail")})
        elif kind == "agent_answered":
            ViewerAgentPort._hear(listener, "answered", {"allowed": event.get("allowed") is True})

    def _end(self, reason: str) -> None:
        with self._lock:
            self._stopped = reason
            calls = list(self._calls.values())
            self._calls.clear()
            slots = list(self._lists.values())
        for listener in calls:
            ViewerAgentPort._hear(listener, "stopped", {"reason": reason})
        for slot in slots:
            if slot[1] is None:
                slot[1] = ("failed", f"The Cellucid viewer stopped ({reason}).")
                slot[0].set()
        self._closed.set()

    @property
    def stopped(self) -> str | None:
        return self._stopped
